version 1.0

task perbase_to_pon_vcf {
    input {
        Array[File] perbase_tsvs
        Array[String] samples
        File reference_fasta
        File reference_fasta_index

        String pon_vcf_name = "pon.vcf.gz"
        Int min_samples_covered = 2
        Int min_samples_alt = 2
        Array[String] chroms = [
            "chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8",
            "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16",
            "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY"
        ]

        String docker_image = "jiminpark/perbase-pon-vcf:1.3"
        Int threads = 2
        Int memSizeGB = 8
        Int diskSizeGB = round(size(perbase_tsvs, "GB")) + 20
    }

    command <<<
        set -o pipefail
        set -e
        set -u
        set -o xtrace

        # pysam needs the .fai next to the FASTA; localized files may land in different dirs
        ln -s ~{reference_fasta} reference.fa
        ln -s ~{reference_fasta_index} reference.fa.fai

        python3 /usr/local/bin/perbase_to_pon_vcf.py \
            --tsv ~{sep=' ' perbase_tsvs} \
            --samples ~{sep=' ' samples} \
            --output ~{pon_vcf_name} \
            --reference reference.fa \
            --min-samples-covered ~{min_samples_covered} \
            --min-samples-alt ~{min_samples_alt} \
            --chroms ~{sep=' ' chroms}
    >>>

    output {
        File pon_vcf = "~{pon_vcf_name}"
        File pon_vcf_index = "~{pon_vcf_name}.tbi"
    }

    runtime {
        preemptible: 2
        docker: docker_image
        cpu: threads
        memory: memSizeGB + " GB"
        disks: "local-disk " + diskSizeGB + " SSD"
    }
}
