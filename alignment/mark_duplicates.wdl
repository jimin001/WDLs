version 1.0

# Mark PCR/optical duplicates (flag 0x400) in a coordinate-sorted Illumina BAM, e.g. the
# output of bwa_mem2_ora_chunked.wdl, which does not mark duplicates. Reads are only
# flagged, not removed.
#
# region unset: mark the whole BAM -> <sample>.markdup.bam + index + stats
# region set (e.g. "chr22"): only mark reads in that region, to cheaply estimate the
#   duplicate rate before committing to a whole-genome run. Needs bam_index. Pairs whose
#   mate falls outside the region are counted as single-end reads, so the rate is
#   approximate; a whole chromosome keeps those to a minimum. The full BAM is still
#   localized, so this saves compute, not disk/localization.
workflow markDuplicates {
    input {
        String sample_name
        File bam
        File? bam_index
        String? region
    }

    call markDuplicatesSamtools {
        input:
            sample_name = sample_name,
            bam = bam,
            bam_index = bam_index,
            region = region
    }

    output {
        File markdup_bam = markDuplicatesSamtools.markdup_bam
        File markdup_bam_idx = markDuplicatesSamtools.markdup_bam_idx
        File markdup_stats = markDuplicatesSamtools.markdup_stats
        Float duplicate_rate = markDuplicatesSamtools.duplicate_rate
        Float optical_duplicate_rate = markDuplicatesSamtools.optical_duplicate_rate
    }

}

task markDuplicatesSamtools {
    input {
        # samtools 1.17
        String docker_image = "jiminpark/aligners"

        String sample_name

        File bam
        File? bam_index
        String? region

        # max pixel distance for optical duplicates:
        # 2500 for patterned flow cells (NovaSeq, HiSeq X/4000), 100 for unpatterned (HiSeq 2500, NextSeq)
        Int optical_distance = 2500

        Int threads = 16
        Int memSizeGB = 48
        # a whole-genome run is several hours and restarts from scratch if preempted
        Int preemptible = if defined(region) then 2 else 0
    }
    Int bam_size = ceil(size(bam, "GB"))
    # whole BAM: input + collate temp files + sort temp files + output
    Int diskSizeGB = if defined(region) then bam_size + 50 else ceil(bam_size * 4.5) + 50

    # "chr22:1-1000000" -> "chr22_1-1000000"
    String out_prefix = if defined(region) then sample_name + "." + sub(select_first([region, ""]), ":", "_") else sample_name
    String out_bam = out_prefix + ".markdup.bam"

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

        # markdup only finds duplicates among reads it sees together; a coordinate sorted bam is expected
        if ! samtools view -H ~{bam} | grep -q "^@HD.*SO:coordinate"
        then
            echo "ERROR: bam header is not SO:coordinate" >&2
            exit 1
        fi

        if [[ "~{region}" != "" ]]
        then
            if [[ "~{bam_index}" == "" ]]
            then
                echo "ERROR: bam_index is required when region is set" >&2
                exit 1
            fi
            # index may be localized to a different directory than the bam
            ln -s ~{bam} input.bam
            ln -s ~{bam_index} input.bam.bai
            # primary alignments of paired reads in the region only
            samtools view -@ ~{threads} -u -f 1 -F 0x900 input.bam ~{region} > input.region.bam
            INPUT_BAM=input.region.bam
        else
            INPUT_BAM=~{bam}
        fi

        # collate groups mates together so fixmate -m can add the MC/ms tags markdup needs,
        # then back to coordinate order for markdup
        samtools collate -@ ~{threads} -O -u ${INPUT_BAM} collate_tmp \
        | samtools fixmate -@ ~{threads} -m -u - - \
        | samtools sort -@ ~{threads} -m 2G -T sort_tmp -u - \
        | samtools markdup -@ ~{threads} -S -d ~{optical_distance} -f ~{out_prefix}.markdup_stats.txt - ~{out_bam}

        samtools index -@ ~{threads} ~{out_bam}

        # duplicate rate = duplicate primary reads / reads examined (PAIRED and DUPLICATE PAIR count reads, not pairs)
        awk -F': ' '
            $1 == "EXAMINED" {examined = $2}
            $1 == "DUPLICATE PRIMARY TOTAL" {dup = $2}
            $1 == "DUPLICATE PAIR OPTICAL" {opt += $2}
            $1 == "DUPLICATE SINGLE OPTICAL" {opt += $2}
            END {
                if (examined == 0) {print "ERROR: no reads examined" > "/dev/stderr"; exit 1}
                printf "%.6f\n", dup / examined > "duplicate_rate.txt"
                printf "%.6f\n", opt / examined > "optical_duplicate_rate.txt"
            }' ~{out_prefix}.markdup_stats.txt

        rm -f input.region.bam
    >>>

    output {
        File markdup_bam = "working/~{out_bam}"
        File markdup_bam_idx = "working/~{out_bam}.bai"
        File markdup_stats = "working/~{out_prefix}.markdup_stats.txt"
        Float duplicate_rate = read_float("working/duplicate_rate.txt")
        Float optical_duplicate_rate = read_float("working/optical_duplicate_rate.txt")
    }

    runtime {
        memory: memSizeGB + " GB"
        cpu: threads
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: docker_image
        preemptible: preemptible
    }
}
