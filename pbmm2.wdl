version 1.0

workflow alignHiFiReads {
    call pbmm2

    output {
        File aligned_bam = pbmm2.aligned_bam
        File aligned_bam_idx = pbmm2.aligned_bam_idx
    }
}

task pbmm2 {
    input {
        String docker_image = "quay.io/biocontainers/pbmm2:1.13.1--h9ee0642_0"

        File ubam1
        File ubam2

        File reference
        String sample_name

        String? additional_args

        Int threads = 64
        Int memSizeGB = 64
    }

    Int file_size = ceil(size(ubam1, "GB"))
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

        # check if length of "additional_args" is zero

        if [[ "~{additional_args}" == "" ]]
        then
            ADDITIONAL_ARGS=""
        else
            ADDITIONAL_ARGS="~{additional_args}"
        fi

        mkdir working
        cd working

        echo ~{ubam1} > myfiles.fofn
        echo ~{ubam2} >> myfiles.fofn

        pbmm2 align ~{reference} myfiles.fofn ~{sample_name}_HiFi.GRCh38.sorted.bam \
        --sort \
        --min-length 50 \
        --sample ~{sample_name} \
        --preset HiFi \
        --num-threads ~{threads} \
        -J 4 \
        --log-level INFO \
        ${ADDITIONAL_ARGS}
    >>>

    output {
        File aligned_bam = "working/~{sample_name}_HiFi.GRCh38.sorted.bam"
        File aligned_bam_idx = "working/~{sample_name}_HiFi.GRCh38.sorted.bam.bai"
    }

    runtime {
        memory: memSizeGB + " GB"
        cpu: threads
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: docker_image
    }
}
