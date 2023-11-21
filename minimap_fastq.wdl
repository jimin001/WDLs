version 1.0

workflow alignLongReads {
    input {
        Array[File] fastq_files
        String sample_name
        String data_type

    }

    scatter (fastq in fastq_files) {
        call minimap2 {
            input:
                fastq = fastq
        }
    }

    call mergeBam {
        input:
            sample_name = sample_name,
            data_type = data_type,
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

        File fastq
        File reference
        String temp_prefix = basename(fastq, ".fastq.gz")
        String? additional_args

        Int minimap_threads
        Int samtools_threads

        Int memSizeGB = 64
        
    }

    Int file_size = ceil(size(fastq, "GB"))
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

        minimap2 -ax map-ont \
        -t ~{minimap_threads} \
        ${ADDITIONAL_ARGS} \
        ~{reference} \
        ~{fastq} \
        | samtools sort -@ ~{samtools_threads} -m 4G > ~{temp_prefix}_GRCh38.sorted.bam

    >>>

    output {
        Array[File] aligned_bam = glob("working/*.bam")
    }

    runtime {
        memory: memSizeGB + " GB"
        cpu: minimap_threads
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: docker_image
    }
}

task mergeBam {
    input {
        String docker_image = "jiminpark/aligners"

        Array[File] aligned_bams
        String sample_name
        String data_type

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

        samtools merge -o "~{sample_name}.~{data_type}.GRCh38.sorted.bam" ~{sep=" " aligned_bams}
        samtools index -@ ~{threads} "~{sample_name}.~{data_type}.GRCh38.sorted.bam"

    >>>

    output {
        File merged_bam = "~{sample_name}.~{data_type}.GRCh38.sorted.bam"
        File merged_bam_idx = "~{sample_name}.~{data_type}.GRCh38.sorted.bam.bai"
    }

    runtime {
        memory: memSizeGB + " GB"
        cpu: threads
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: docker_image
    }
}
