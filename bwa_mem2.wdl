version 1.0

workflow alignShortReads {

    call bwaMem2

    output {
        File aligned_bam = bwaMem2.aligned_bam
        File aligned_bam_idx = bwaMem2.aligned_bam_idx
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
    }
    Int file_size = ceil(size(fastq1, "GB"))
    Int diskSizeGB = 3 * file_size

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
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: docker_image
    }
}


