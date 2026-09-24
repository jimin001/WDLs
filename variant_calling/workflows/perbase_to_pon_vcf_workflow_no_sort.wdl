version 1.0

# Same as perbase_to_pon_vcf_workflow.wdl but skips the sort_perbase_tsv
# check/sort step entirely -- perbase_tsvs are passed straight into
# perbase_to_pon_vcf, which still assumes they're sorted per-chromosome in
# the given --chroms order and will raise if that assumption is violated.

import "../tasks/perbase_to_pon_vcf_task.wdl" as perbase_to_pon_vcf_task

workflow PerbaseToPonVcf {
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
        Int diskSizeGB = 0
    }

    call perbase_to_pon_vcf_task.perbase_to_pon_vcf as perbase_to_pon_vcf {
        input:
            perbase_tsvs = perbase_tsvs,
            samples = samples,
            reference_fasta = reference_fasta,
            reference_fasta_index = reference_fasta_index,
            pon_vcf_name = pon_vcf_name,
            min_samples_covered = min_samples_covered,
            min_samples_alt = min_samples_alt,
            chroms = chroms,
            docker_image = docker_image,
            threads = threads,
            memSizeGB = memSizeGB,
            diskSizeGB = if diskSizeGB > 0 then diskSizeGB else round(size(perbase_tsvs, "GB")) + 20
    }

    output {
        File pon_vcf = perbase_to_pon_vcf.pon_vcf
        File pon_vcf_index = perbase_to_pon_vcf.pon_vcf_index
    }
}
