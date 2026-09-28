version 1.0

task perbase_to_pon_vcf {
    input {
        Array[File] perbase_tsvs
        Array[String] samples
        # Optional tabix indexes (index_perbase_tsv), 1:1 with perbase_tsvs.
        # When given, each TSV is read by fetching only the --chroms contigs.
        Array[File] perbase_tsv_indexes = []
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
        # Contigs declared in the VCF header (default: chroms). Set to the full
        # contig list when chroms is one shard of a per-chromosome split.
        Array[String] header_chroms = []

        String docker_image = "jiminpark/perbase-pon-vcf:1.4"
        Int threads = 2
        Int memSizeGB = 8
        Int diskSizeGB = round(size(perbase_tsvs, "GB") + size(perbase_tsv_indexes, "GB")) + 20
    }

    command <<<
        set -o pipefail
        set -e
        set -u
        set -o xtrace

        # pysam needs the .fai next to the FASTA; localized files may land in different dirs
        ln -s ~{reference_fasta} reference.fa
        ln -s ~{reference_fasta_index} reference.fa.fai

        tsvs=(~{sep=' ' perbase_tsvs})
        indexes=(~{sep=' ' perbase_tsv_indexes})
        if [ ${#indexes[@]} -gt 0 ]; then
            if [ ${#indexes[@]} -ne ${#tsvs[@]} ]; then
                echo "ERROR: perbase_tsv_indexes must pair 1:1 with perbase_tsvs" >&2
                exit 1
            fi
            # the script finds each TSV's index as <tsv>.tbi, and localized
            # files may land in different dirs, so pair them up side by side
            for i in "${!tsvs[@]}"; do
                ln -s "${tsvs[i]}" "input_${i}.tsv.gz"
                ln -s "${indexes[i]}" "input_${i}.tsv.gz.tbi"
                tsvs[i]="input_${i}.tsv.gz"
            done
        fi

        header_chroms=(~{sep=' ' header_chroms})
        header_args=()
        if [ ${#header_chroms[@]} -gt 0 ]; then
            header_args=(--header-chroms "${header_chroms[@]}")
        fi

        python3 /usr/local/bin/perbase_to_pon_vcf.py \
            --tsv "${tsvs[@]}" \
            --samples ~{sep=' ' samples} \
            --output ~{pon_vcf_name} \
            --reference reference.fa \
            --min-samples-covered ~{min_samples_covered} \
            --min-samples-alt ~{min_samples_alt} \
            --chroms ~{sep=' ' chroms} \
            "${header_args[@]}"
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
