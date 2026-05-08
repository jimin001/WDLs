version 1.0

task check_pacbio_signal {
    input {
        File ont_vcf
        File ont_vcf_idx
        File pacbio_bam
        File pacbio_bai
        File script

        String sample

        String mode = "snv"
        Int min_mapq = 20
        Int min_baseq = 10
        Int min_alt_count = 1
        Float min_alt_freq = 0.0

        String docker_image = "jiminpark/deepsomatic_postprocess:v3"
        Int threads = 8
        Int memSizeGB = 16
        Int diskSizeGB = round(size(pacbio_bam, "G")) * 3
    }

    command <<<
        set -o pipefail
        set -e
        set -u
        set -o xtrace

        python3 ~{script} \
            -vcf ~{ont_vcf} \
            -bam ~{pacbio_bam} \
            -o ~{sample}_pacbio_signal_output.tsv \
            --mode ~{mode} \
            -q ~{min_mapq} \
            -Q ~{min_baseq} \
            --min-alt-count ~{min_alt_count} \
            --min-alt-freq ~{min_alt_freq} \
            -t ~{threads}
    >>>

    output {
        File signal_tsv = "~{sample}_pacbio_signal_output.tsv"
    }

    runtime {
        docker: docker_image
        cpu: threads
        memory: memSizeGB + " GB"
        disks: "local-disk " + diskSizeGB + " SSD"
    }
}
