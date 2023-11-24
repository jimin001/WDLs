version 1.0

workflow alignONTReads {
    input {
        Array[File] ubam_files
        String sample_name

    }

    scatter (ubam in ubam_files) {
        call minimap2 {
            input:
                ubam = ubam
        }
    }

    call mergeBam {
        input:
            sample_name = sample_name,
            aligned_bams = flatten(minimap2.aligned_bam)

    }

    output {
        File merged_bam = mergeBam.merged_bam
        File merged_bam_idx = mergeBam.merged_bam_idx
    }
}


task minimap2 {
    input {
        String docker_image = "jiminpark/aligners"

        File ubam
        File reference
        String temp_prefix = basename(ubam, ".bam")
        String? additional_args

        Int threads = 64
        Int memSizeGB = 64
        
    }

    Int file_size = ceil(size(ubam, "GB"))
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

        samtools fastq -T MM,ML ~{ubam} | minimap2 -ax map-ont ~{reference} - \
        -t ~{threads} \
        ${ADDITIONAL_ARGS} \
        | samtools sort -@4 -m 4G > ~{temp_prefix}_GRCh38.sorted.bam

    >>>

    output {
        Array[File] aligned_bam = glob("working/*.bam")
    }

    runtime {
        memory: memSizeGB + " GB"
        cpu: threads
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: docker_image
    }
}

task mergeBam {
    input {
        String docker_image = "jiminpark/aligners"

        Array[File] aligned_bams
        String sample_name

        Int threads = 10
        Int memSizeGB = 64

    }

    Int file_size = ceil(size(aligned_bams[0], "GB"))
    Int diskSizeGB = 4 * file_size


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

        samtools merge -o "~{sample_name}_ONT.GRCh38.sorted.bam" ~{sep=" " aligned_bams}
        samtools index -@ ~{threads} "~{sample_name}_ONT.GRCh38.sorted.bam"

    >>>

    output {
        File merged_bam = "~{sample_name}_ONT.GRCh38.sorted.bam"
        File merged_bam_idx = "~{sample_name}_ONT.GRCh38.sorted.bam.bai"
    }

    runtime {
        memory: memSizeGB + " GB"
        cpu: threads
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: docker_image
    }
}
