version 1.0

task perbase {
    input {
        File bam
        File bam_idx
        String sample

        File? bed_file
        Int min_base_quality_score = 20
        Int min_mapq = 10

        String docker_image = "jiminpark/perbase:1.4.0"
        Int threads = 8
        Int memSizeGB = 32
        Int diskSizeGB = 3 * round(size(bam, "G")) + 50
    }

    command <<<
        set -o pipefail
        set -e
        set -u
        set -o xtrace

        ln -s ~{bam} ~{sample}.bam
        ln -s ~{bam_idx} ~{sample}.bam.bai

        perbase base-depth ~{sample}.bam \
            -o ~{sample}_perbase_output.tsv \
            -t ~{threads} \
            --min-base-quality-score ~{min_base_quality_score} \
            --min-mapq ~{min_mapq} \
            ~{"-b " + bed_file}
    >>>

    output {
        File perbase_output = "~{sample}_perbase_output.tsv"
    }

    runtime {
        preemptible: 2
        docker: docker_image
        cpu: threads
        memory: memSizeGB + " GB"
        disks: "local-disk " + diskSizeGB + " SSD"
    }
}
