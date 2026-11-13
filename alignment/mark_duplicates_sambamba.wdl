version 1.0

# Mark duplicates (flag 0x400) in an already aligned, merged, coordinate-sorted Illumina BAM,
# e.g. the output of bwa_mem2_ora_chunked.wdl. Reads are only flagged, not removed.
#
# Cheaper alternative to mark_duplicates.wdl for whole-genome runs: sambamba markdup reads the
# coordinate-sorted BAM directly, so the whole-genome samtools collate + fixmate + sort that
# samtools markdup needs is skipped (less runtime, about half the disk).
#
# Differences from mark_duplicates.wdl:
#   - no optical vs PCR duplicate split (use mark_duplicates.wdl region mode on e.g. chr22 for that)
#   - supplementary alignments of duplicate reads are not marked (bwa-mem2 -M writes split
#     hits as secondary, so bwa_mem2_ora_chunked.wdl BAMs have no supplementary alignments)
#   - no region mode
#
# markDuplicatesSambamba: sambamba markdup + index (non-preemptible, the long step)
# duplicateStats:         samtools flagstat -> duplicate rate, and optionally mean depth with
#                         and without duplicates (preemptible, separate image with samtools)
workflow markDuplicatesSambambaWorkflow {
    input {
        String sample_name
        File bam
        Boolean calculate_coverage = true
    }

    call markDuplicatesSambamba {
        input:
            sample_name = sample_name,
            bam = bam
    }

    call duplicateStats {
        input:
            sample_name = sample_name,
            bam = markDuplicatesSambamba.markdup_bam,
            bam_index = markDuplicatesSambamba.markdup_bam_idx,
            calculate_coverage = calculate_coverage
    }

    output {
        File markdup_bam = markDuplicatesSambamba.markdup_bam
        File markdup_bam_idx = markDuplicatesSambamba.markdup_bam_idx
        File markdup_log = markDuplicatesSambamba.markdup_log
        File flagstat = duplicateStats.flagstat
        Float duplicate_rate = duplicateStats.duplicate_rate
        File? coverage_report = duplicateStats.coverage_report
        Float? mean_coverage_with_duplicates = duplicateStats.mean_coverage_with_duplicates
        Float? mean_coverage_without_duplicates = duplicateStats.mean_coverage_without_duplicates
    }
}

task markDuplicatesSambamba {
    input {
        String docker_image = "quay.io/biocontainers/sambamba:1.0.1--he614052_4"

        String sample_name
        File bam

        # reads held while waiting for their mate; should be > coverage * insert size.
        # Larger tables mean fewer temp files on deep WGS at the cost of memory
        Int hash_table_size = 1048576
        Int overflow_list_size = 1000000
        Int sort_buffer_size_mb = 8192
        Int io_buffer_size_mb = 512

        Int threads = 16
        Int memSizeGB = 32
        # one long pass that restarts from scratch if preempted
        Int preemptible = 0
    }
    Int bam_size = ceil(size(bam, "GB"))
    # input + output + temp files (a fraction of the BAM)
    Int diskSizeGB = ceil(bam_size * 2.5) + 50

    String out_bam = sample_name + ".markdup.bam"

    command <<<
        # Set the exit code of a pipeline to that of the rightmost command
        # to exit with a non-zero status, or zero if all commands of the pipeline exit
        set -o pipefail
        # cause a bash script to exit immediately when a command fails
        set -e
        # cause the bash shell to treat unset variables as an error and exit immediately
        set -u
        # echo each line of the script to stdout so we can see what is happening
        # to turn off echo do 'set +o xtrace'
        set -o xtrace

        mkdir working tmp
        cd working

        # sambamba markdup pairs reads by streaming through coordinate order
        if ! sambamba view -H ~{bam} | grep -q "^@HD.*SO:coordinate"
        then
            echo "ERROR: bam header is not SO:coordinate" >&2
            exit 1
        fi

        # many overflow temp files can be open at once on deep data
        ulimit -n $(ulimit -Hn) || true

        sambamba markdup \
            -t ~{threads} \
            --tmpdir=../tmp \
            --hash-table-size=~{hash_table_size} \
            --overflow-list-size=~{overflow_list_size} \
            --sort-buffer-size=~{sort_buffer_size_mb} \
            --io-buffer-size=~{io_buffer_size_mb} \
            ~{bam} ~{out_bam} \
            2> ~{sample_name}.markdup.log \
        || { cat ~{sample_name}.markdup.log >&2; exit 1; }
        cat ~{sample_name}.markdup.log >&2

        # markdup writes an index alongside its output; build one if this version did not
        if [[ ! -s ~{out_bam}.bai ]]
        then
            sambamba index -t ~{threads} ~{out_bam} ~{out_bam}.bai
        fi

        rm -rf ../tmp
    >>>

    output {
        File markdup_bam = "working/~{out_bam}"
        File markdup_bam_idx = "working/~{out_bam}.bai"
        File markdup_log = "working/~{sample_name}.markdup.log"
    }

    runtime {
        memory: memSizeGB + " GB"
        cpu: threads
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: docker_image
        preemptible: preemptible
    }
}

task duplicateStats {
    input {
        # samtools 1.17
        String docker_image = "jiminpark/aligners"

        String sample_name
        File bam
        File bam_index
        Boolean calculate_coverage = true

        Int threads = 16
        Int memSizeGB = 16
        # only reads the marked BAM, so a preemption just repeats this cheap step
        Int preemptible = 2

        # never set: WDL 1.0 has no None literal for the optional coverage outputs
        File? NONE_FILE
        Float? NONE_FLOAT
    }
    Int diskSizeGB = ceil(size(bam, "GB")) + 30

    command <<<
        # Set the exit code of a pipeline to that of the rightmost command
        # to exit with a non-zero status, or zero if all commands of the pipeline exit
        set -o pipefail
        # cause a bash script to exit immediately when a command fails
        set -e
        # cause the bash shell to treat unset variables as an error and exit immediately
        set -u
        # echo each line of the script to stdout so we can see what is happening
        # to turn off echo do 'set +o xtrace'
        set -o xtrace

        mkdir working
        cd working

        # index may be localized to a different directory than the bam
        ln -s ~{bam} input.bam
        ln -s ~{bam_index} input.bam.bai

        samtools flagstat -@ ~{threads} -O tsv input.bam > ~{sample_name}.flagstat.tsv

        # duplicate rate = primary duplicates / primary mapped reads (QC-passed column),
        # comparable to EXAMINED-based rate from samtools markdup in mark_duplicates.wdl
        awk -F'\t' '
            $3 == "primary mapped" {mapped = $1}
            $3 == "primary duplicates" {dup = $1}
            END {
                if (mapped == 0) {print "ERROR: no primary mapped reads" > "/dev/stderr"; exit 1}
                printf "%.6f\n", dup / mapped > "duplicate_rate.txt"
            }' ~{sample_name}.flagstat.tsv

        if [[ "~{calculate_coverage}" == "true" ]]
        then
            samtools idxstats input.bam | cut -f 1 | grep -E '^(chr)?([0-9]+|X|Y)$' > contigs.txt

            # samtools coverage is single threaded, so run each contig (with and without duplicates) in parallel.
            # --ff defaults to UNMAP,SECONDARY,QCFAIL,DUP; dropping DUP counts the duplicates
            mkdir coverage_tmp
            while read -r CONTIG
            do
                echo "UNMAP,SECONDARY,QCFAIL ${CONTIG} coverage_tmp/with_dups.${CONTIG}.tsv"
                echo "UNMAP,SECONDARY,QCFAIL,DUP ${CONTIG} coverage_tmp/without_dups.${CONTIG}.tsv"
            done < contigs.txt \
            | xargs -P ~{threads} -n 3 sh -c 'samtools coverage -H -d 0 --ff "$0" -r "$1" input.bam > "$2"'

            echo -e "contig\tstart\tend\tlength\tmeandepth_with_duplicates\tmeandepth_without_duplicates" > ~{sample_name}.coverage.tsv
            while read -r CONTIG
            do
                paste coverage_tmp/with_dups.${CONTIG}.tsv coverage_tmp/without_dups.${CONTIG}.tsv \
                | awk -v OFS='\t' '{print $1, $2, $3, $3 - $2 + 1, $7, $16}'
            done < contigs.txt >> ~{sample_name}.coverage.tsv

            # length-weighted mean depth across contigs
            awk -F'\t' 'NR > 1 {len += $4; with += $4 * $5; without += $4 * $6}
                END {
                    printf "%.4f\n", with / len > "mean_coverage_with_duplicates.txt"
                    printf "%.4f\n", without / len > "mean_coverage_without_duplicates.txt"
                }' ~{sample_name}.coverage.tsv

            rm -rf coverage_tmp
        fi
    >>>

    output {
        File flagstat = "working/~{sample_name}.flagstat.tsv"
        Float duplicate_rate = read_float("working/duplicate_rate.txt")
        File? coverage_report = if calculate_coverage then "working/~{sample_name}.coverage.tsv" else NONE_FILE
        Float? mean_coverage_with_duplicates = if calculate_coverage then read_float("working/mean_coverage_with_duplicates.txt") else NONE_FLOAT
        Float? mean_coverage_without_duplicates = if calculate_coverage then read_float("working/mean_coverage_without_duplicates.txt") else NONE_FLOAT
    }

    runtime {
        memory: memSizeGB + " GB"
        cpu: threads
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: docker_image
        preemptible: preemptible
    }
}
