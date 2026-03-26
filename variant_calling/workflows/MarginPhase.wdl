version 1.0

workflow MarginPhase {
    input {
        File bam
        File bam_index
        File reference
        File reference_index
        File vcf
        File vcf_index
        String output_prefix
        Int threads = 64
    }

    call marginPhase {
        input:
            bam            = bam,
            bam_index      = bam_index,
            reference      = reference,
            reference_index = reference_index,
            vcf            = vcf,
            vcf_index      = vcf_index,
            output_prefix  = output_prefix,
            threads        = threads
    }

    call indexBam {
        input:
            bam     = marginPhase.phased_bam,
            threads = threads
    }

    output {
        File phased_bam_index           = indexBam.bam_index
        Array[File] phased_output_files = marginPhase.phased_output_files
    }
}


task marginPhase {
    input {
        String docker_image = "meredith705/card_harmonize_vcf:0.2"

        File bam
        File bam_index
        File reference
        File reference_index
        File vcf
        File vcf_index
        String output_prefix
        String params_json = "/opt/margin/params/phase/allParams.phase_vcf.ont.json"

        Int threads     = 64
        Int memSizeGB   = 200
        Int diskSizeGB  = 512
    }

    command <<<
        set -o pipefail
        set -e
        set -u
        set -o xtrace

        margin phase \
            ~{bam} \
            ~{reference} \
            ~{vcf} \
            ~{params_json} \
            -t ~{threads} \
            -o ~{output_prefix}

    >>>

    output {
        File phased_bam                 = "~{output_prefix}.phased.bam"
        Array[File] phased_output_files = glob("~{output_prefix}*")
    }

    runtime {
        docker: docker_image
        memory: memSizeGB + " GB"
        cpu: threads
        disks: "local-disk " + diskSizeGB + " SSD"
    }
}


task indexBam {
    input {
        String docker_image = "jiminpark/deepsomatic_postprocess:v3"

        File bam
        Int threads   = 30
        Int memSizeGB = 16
        Int diskSizeGB = 256
    }

    command <<<
        set -o pipefail
        set -e
        set -u
        set -o xtrace

        samtools index -@ ~{threads} ~{bam}
    >>>

    output {
        File bam_index = "~{bam}.bai"
    }

    runtime {
        docker: docker_image
        memory: memSizeGB + " GB"
        cpu: threads
        disks: "local-disk " + diskSizeGB + " SSD"
    }
}
