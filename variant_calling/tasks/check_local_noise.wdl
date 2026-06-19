version 1.0

task check_local_noise {
    input {
        File bam
        File bai
        File vcf
        File vcf_idx
        File ref
        File ref_fai
        File script

        String sample

        Boolean tag_vcf = false
        Int window = 10
        Int min_mapq = 20

        String docker_image = "jiminpark/deepsomatic_postprocess:v3"
        Int memSizeGB = 16
        Int diskSizeGB = round(size(bam, "G")) * 2 + round(size(ref, "G")) + 10
    }

    command <<<
        set -o pipefail
        set -e
        set -u
        set -o xtrace

        python3 ~{script} \
            ~{bam} \
            ~{vcf} \
            ~{ref} \
            --window ~{window} \
            --min-mapq ~{min_mapq} \
            -o ~{sample}.local_cleanliness.tsv \
            ~{if tag_vcf then "--vcf-out " + sample + ".local_cleanliness.vcf.gz" else ""}

        if ~{tag_vcf}; then
            tabix -p vcf ~{sample}.local_cleanliness.vcf.gz
        fi
    >>>

    output {
        File        cleanliness_tsv = "~{sample}.local_cleanliness.tsv"
        Array[File] tagged_vcf      = glob("~{sample}.local_cleanliness.vcf.gz")
        Array[File] tagged_vcf_idx  = glob("~{sample}.local_cleanliness.vcf.gz.tbi")
    }

    runtime {
        docker: docker_image
        memory: memSizeGB + " GB"
        disks: "local-disk " + diskSizeGB + " SSD"
    }
}
