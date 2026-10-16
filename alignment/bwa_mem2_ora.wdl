version 1.0

workflow alignShortReadsOra {

    input {
        String sample_name

        File fastq1_ora
        File fastq2_ora

        File reference

        String? additional_args

        # override the computed disk size for bwaMem2
        Int? bwa_diskSizeGB
    }

    call decompressOra as decompressOra1 {
        input:
            ora_fastq = fastq1_ora
    }

    call decompressOra as decompressOra2 {
        input:
            ora_fastq = fastq2_ora
    }

    call bwaMem2 {
        input:
            sample_name = sample_name,
            fastq1 = decompressOra1.fastq_gz,
            fastq2 = decompressOra2.fastq_gz,
            reference = reference,
            additional_args = additional_args,
            diskSizeGB = bwa_diskSizeGB
    }

    output {
        File aligned_bam = bwaMem2.aligned_bam
        File aligned_bam_idx = bwaMem2.aligned_bam_idx
    }

}

task decompressOra {
    input {
        # built from Dockerfile.orad (orad on PATH, ORA_REF_PATH set to the bundled oradata/)
        String docker_image = "jiminpark/orad:2.7.0"

        File ora_fastq

        Int threads = 16
        Int memSizeGB = 32
        Int preemptible = 2
    }
    # ORA is ~5x smaller than fastq.gz; leave room for the gz output
    Int ora_size = ceil(size(ora_fastq, "GB"))
    Int diskSizeGB = 8 * ora_size + 20

    # e.g. sample_R1_001.fastq.ora -> sample_R1_001.fastq.gz
    String out_name = basename(ora_fastq, ".fastq.ora") + ".fastq.gz"

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

        # decompress ora to fastq.gz
        orad \
            -t ~{threads} \
            -o working/~{out_name} \
            ~{ora_fastq}

        # quick sanity check: ora file stats + first two reads of the output
        orad -i ~{ora_fastq} > working/~{out_name}.summary.txt
        echo "--- first 2 reads ---" >> working/~{out_name}.summary.txt
        { zcat working/~{out_name} || true; } | head -n 8 >> working/~{out_name}.summary.txt
    >>>

    output {
        File fastq_gz = "working/~{out_name}"
        File summary = "working/~{out_name}.summary.txt"
    }

    runtime {
        memory: memSizeGB + " GB"
        cpu: threads
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: docker_image
        preemptible: preemptible
    }
}

task bwaMem2 {
    input {
        String docker_image = "jiminpark/aligners"

        String sample_name

        File fastq1
        File fastq2

        File reference

        String? additional_args

        Int threads = 64
        Int memSizeGB = 128

        Int? diskSizeGB
    }
    # both fastqs + sorted BAM + samtools sort temp files + reference/index
    Int fastq_size = ceil(size(fastq1, "GB") + size(fastq2, "GB"))
    Int ref_size = ceil(size(reference, "GB"))
    Int computedDiskSizeGB = ceil(2.5 * fastq_size) + 10 * ref_size + 100

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

        # check if length of "additional_args" is zero

        if [[ "~{additional_args}" == "" ]]
        then
            ADDITIONAL_ARGS=""
        else
            ADDITIONAL_ARGS="~{additional_args}"
        fi

        # index reference
        bwa-mem2 index ~{reference}

        # align with bwa-mem2 and sort with samtools
        bwa-mem2 mem -M ${ADDITIONAL_ARGS} -t ~{threads} ~{reference} ~{fastq1} ~{fastq2} \
        | samtools sort -@4 -m 4G > ~{sample_name}_Illumina.GRCh38.sorted.bam

        # index bam
        samtools index -@ ~{threads} ~{sample_name}_Illumina.GRCh38.sorted.bam
    >>>

    output {
        File aligned_bam = "working/~{sample_name}_Illumina.GRCh38.sorted.bam"
        File aligned_bam_idx = "working/~{sample_name}_Illumina.GRCh38.sorted.bam.bai"
    }

    runtime {
        memory: memSizeGB + " GB"
        cpu: threads
        disks: "local-disk " + select_first([diskSizeGB, computedDiskSizeGB]) + " SSD"
        docker: docker_image
    }
}
