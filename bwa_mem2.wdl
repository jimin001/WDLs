version 1.0

workflow alignShortReads {
    input {
        String sample_name
        String data_type
    }

    call bwaMem2 {
        input:
            sample_name = sample_name,
            data_type = data_type
    }

    call samtools {
        input:
            sample_name = sample_name,
            data_type = data_type,
            aligned_sam = bwaMem2.aligned_sam
    }

    output {
        File aligned_bam = samtools.aligned_bam
        File aligned_bam_idx = samtools.aligned_bam_idx
    }
    
}

task bwaMem2 {
    input {
        String docker_image = "jiminpark/aligners"

        String sample_name
        String data_type

        File fastq1
        File fastq2

        File reference
        File reference_idx

        String? additional_args

        Int threads = 64
        Int memSizeGB = 128
        Int diskSizeGB = 200
    }

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

        bwa-mem2 index ~{reference}
        bwa-mem2 mem -M ${ADDITIONAL_ARGS} -t ~{threads} ~{reference} ~{fastq1} ~{fastq2} > "~{sample_name}_~{data_type}.sam"
    >>>

    output {
        File aligned_sam = "working/~{sample_name}_~{data_type}.sam"
    }

    runtime {
        memory: memSizeGB + " GB"
        cpu: threads
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: docker_image
    }
}

task samtools {
    input {
        String docker_image = "jiminpark/aligners"
        String sample_name
        String data_type

        File aligned_sam

        Int threads = 64
        Int memSizeGB = 64
        Int diskSizeGB = 200

    }

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

        samtools view -b -@ ~{threads} ~{aligned_sam} | samtools sort -@ ~{threads} > "~{sample_name}_~{data_type}.sorted.bam"
        samtools index -@ ~{threads} "~{sample_name}_~{data_type}.sorted.bam"
    >>>

    output {
        File aligned_bam = "working/~{sample_name}_~{data_type}.sorted.bam"
        File aligned_bam_idx = "working/~{sample_name}_~{data_type}.sorted.bam.bai"
    }

    runtime {
        memory: memSizeGB + " GB"
        cpu: threads
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: docker_image
    }

}
