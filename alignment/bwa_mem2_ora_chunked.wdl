version 1.0

# Chunked version of bwa_mem2_ora.wdl. A single 64-core alignment runs for many hours
# and was preempted after ~6h (twice), restarting from scratch each time. Here the reads
# are split into chunks that each align in ~1h on a smaller preemptible VM, so a
# preemption only loses one chunk.
#
# splitOra:   ORA -> orad --raw --stdout -> split into fastq.gz chunks (R1/R2 in lockstep)
# alignChunk: scatter, bwa-mem2 + samtools sort per chunk
# mergeBams:  samtools merge the sorted chunk BAMs + index
workflow alignShortReadsOraChunked {

    input {
        String sample_name

        File fastq1_ora
        File fastq2_ora

        # prebuilt bwa-mem2 index (see bwa_mem2_index.wdl); the fasta itself is not needed
        File ref_0123
        File ref_amb
        File ref_ann
        File ref_bwt_2bit_64
        File ref_pac

        # read pairs per chunk; ~1h of bwa-mem2 on 32 threads (tune from alignChunk runtimes)
        Int reads_per_chunk = 40000000

        String? additional_args
    }

    call splitOra {
        input:
            fastq1_ora = fastq1_ora,
            fastq2_ora = fastq2_ora,
            reads_per_chunk = reads_per_chunk
    }

    scatter (chunk in zip(splitOra.r1_chunks, splitOra.r2_chunks)) {
        call alignChunk {
            input:
                fastq1 = chunk.left,
                fastq2 = chunk.right,
                ref_0123 = ref_0123,
                ref_amb = ref_amb,
                ref_ann = ref_ann,
                ref_bwt_2bit_64 = ref_bwt_2bit_64,
                ref_pac = ref_pac,
                additional_args = additional_args
        }
    }

    call mergeBams {
        input:
            sample_name = sample_name,
            chunk_bams = alignChunk.chunk_bam
    }

    output {
        File aligned_bam = mergeBams.aligned_bam
        File aligned_bam_idx = mergeBams.aligned_bam_idx
        File ora_summary = splitOra.ora_summary
    }

}

task splitOra {
    input {
        # built from Dockerfile.bwa_mem2_ora (jiminpark/aligners + orad, ORA_REF_PATH set)
        String docker_image = "jiminpark/bwa_mem2_ora:2.7.0"

        File fastq1_ora
        File fastq2_ora

        Int reads_per_chunk

        # per read (R1 and R2 run in parallel): orad threads + bgzip threads
        Int orad_threads = 4
        Int bgzip_threads = 6
        Int threads = 16
        Int memSizeGB = 32
        Int preemptible = 2
    }
    # ORA is ~5x smaller than fastq.gz; level-1 bgzip chunks are somewhat larger than that
    Int ora_size = ceil(size(fastq1_ora, "GB") + size(fastq2_ora, "GB"))
    Int diskSizeGB = 9 * ora_size + 20

    Int lines_per_chunk = 4 * reads_per_chunk

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

        # quick sanity check: ora file stats
        orad -i ~{fastq1_ora} > ora_summary.txt
        orad -i ~{fastq2_ora} >> ora_summary.txt

        # R1 and R2 are split by the same line count, so chunk N of R1 pairs with chunk N of R2.
        # bgzip level 1: chunks are short-lived intermediates, speed matters more than size
        orad --raw --stdout -t ~{orad_threads} ~{fastq1_ora} \
        | split -l ~{lines_per_chunk} -d -a 4 \
            --filter='bgzip -@ ~{bgzip_threads} -l 1 > ${FILE}.fastq.gz' \
            - R1_chunk_ &
        SPLIT1_PID=$!
        orad --raw --stdout -t ~{orad_threads} ~{fastq2_ora} \
        | split -l ~{lines_per_chunk} -d -a 4 \
            --filter='bgzip -@ ~{bgzip_threads} -l 1 > ${FILE}.fastq.gz' \
            - R2_chunk_ &
        SPLIT2_PID=$!

        # collect both before failing (if one side dies the other may be killed too)
        set +e
        wait ${SPLIT1_PID}; SPLIT1_RC=$?
        wait ${SPLIT2_PID}; SPLIT2_RC=$?
        set -e
        if [[ ${SPLIT1_RC} -ne 0 || ${SPLIT2_RC} -ne 0 ]]
        then
            echo "orad/split failed: R1 exit ${SPLIT1_RC}, R2 exit ${SPLIT2_RC}" >&2
            exit 1
        fi

        N_R1=$(ls R1_chunk_*.fastq.gz | wc -l)
        N_R2=$(ls R2_chunk_*.fastq.gz | wc -l)
        if [[ ${N_R1} -ne ${N_R2} ]]
        then
            echo "R1 and R2 produced different numbers of chunks; read counts differ" >&2
            exit 1
        fi
        echo "n_chunks: ${N_R1}" >> ora_summary.txt
    >>>

    output {
        # glob, not read_lines: Cromwell on Terra only delocalizes outputs it can resolve
        # before reading task files. alignChunk checks the R1/R2 chunk numbers match
        Array[File] r1_chunks = glob("working/R1_chunk_*.fastq.gz")
        Array[File] r2_chunks = glob("working/R2_chunk_*.fastq.gz")
        File ora_summary = "working/ora_summary.txt"
    }

    runtime {
        memory: memSizeGB + " GB"
        cpu: threads
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: docker_image
        preemptible: preemptible
    }
}

task alignChunk {
    input {
        String docker_image = "jiminpark/bwa_mem2_ora:2.7.0"

        File fastq1
        File fastq2

        File ref_0123
        File ref_amb
        File ref_ann
        File ref_bwt_2bit_64
        File ref_pac

        String? additional_args

        Int threads = 32
        Int sort_threads = 4
        Int memSizeGB = 64
        # ~1h chunks are rarely preempted, and a preemption only costs that chunk
        Int preemptible = 3
    }
    Int fastq_size = ceil(size(fastq1, "GB") + size(fastq2, "GB"))
    Int index_size = ceil(size([ref_0123, ref_amb, ref_ann, ref_bwt_2bit_64, ref_pac], "GB"))
    Int diskSizeGB = 4 * fastq_size + index_size + 20

    # R1_chunk_0003.fastq.gz -> chunk_0003
    String chunk_id = sub(basename(fastq1, ".fastq.gz"), "^R1_", "")
    String chunk2_id = sub(basename(fastq2, ".fastq.gz"), "^R2_", "")
    String chunk_bam_name = chunk_id + ".sorted.bam"

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

        # R1/R2 come from two separate globs; never align mismatched chunks
        if [[ "~{chunk_id}" != "~{chunk2_id}" ]]
        then
            echo "R1/R2 chunk mismatch: ~{chunk_id} vs ~{chunk2_id}" >&2
            exit 1
        fi

        mkdir working ref
        cd working

        # check if length of "additional_args" is zero

        if [[ "~{additional_args}" == "" ]]
        then
            ADDITIONAL_ARGS=""
        else
            ADDITIONAL_ARGS="~{additional_args}"
        fi

        # index files may be localized to different directories; put them under one prefix
        ln -s ~{ref_0123} ../ref/ref.fa.0123
        ln -s ~{ref_amb} ../ref/ref.fa.amb
        ln -s ~{ref_ann} ../ref/ref.fa.ann
        ln -s ~{ref_bwt_2bit_64} ../ref/ref.fa.bwt.2bit.64
        ln -s ~{ref_pac} ../ref/ref.fa.pac

        # align with bwa-mem2 and sort with samtools
        bwa-mem2 mem -M ${ADDITIONAL_ARGS} -t ~{threads} ../ref/ref.fa ~{fastq1} ~{fastq2} \
        | samtools sort -@ ~{sort_threads} -m 2G -o ~{chunk_bam_name} -
    >>>

    output {
        File chunk_bam = "working/~{chunk_bam_name}"
    }

    runtime {
        memory: memSizeGB + " GB"
        cpu: threads
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: docker_image
        preemptible: preemptible
    }
}

task mergeBams {
    input {
        String docker_image = "jiminpark/bwa_mem2_ora:2.7.0"

        String sample_name

        Array[File] chunk_bams

        Int threads = 16
        Int memSizeGB = 16
        Int preemptible = 2
    }
    Int bam_size = ceil(size(chunk_bams, "GB"))
    Int diskSizeGB = 2 * bam_size + 20

    String out_bam = "~{sample_name}_Illumina.GRCh38.sorted.bam"

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

        # chunk BAMs are each coordinate-sorted, so this is a single streaming merge
        samtools merge -@ ~{threads} -b ~{write_lines(chunk_bams)} -o ~{out_bam}

        # index bam
        samtools index -@ ~{threads} ~{out_bam}
    >>>

    output {
        File aligned_bam = "working/~{out_bam}"
        File aligned_bam_idx = "working/~{out_bam}.bai"
    }

    runtime {
        memory: memSizeGB + " GB"
        cpu: threads
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: docker_image
        preemptible: preemptible
    }
}
