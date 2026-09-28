version 1.0

# Same output as perbase_to_pon_vcf_workflow_no_sort.wdl, but split by
# chromosome so wall time is bounded by the largest chromosome instead of the
# whole genome:
#   1. index_perbase_tsv: tabix-index each sample's bgzipped perbase TSV, in
#      parallel. This also enforces the sortedness perbase_to_pon_vcf relies
#      on (tabix fails on unsorted input), so no sort step is needed.
#   2. perbase_to_pon_vcf: one shard per chromosome, each fetching only its
#      chromosome from every sample via the index.
#   3. concat_pon_vcf: join the shards in chroms order into one VCF.
#
# perbase_tsvs must be bgzip-compressed (as the Perbase WDL writes them).
# Every shard localizes every sample's full TSV, so each shard needs disk for
# the whole panel even though it reads only its own chromosome.

import "../tasks/index_perbase_tsv_task.wdl" as index_perbase_tsv_task
import "../tasks/perbase_to_pon_vcf_task.wdl" as perbase_to_pon_vcf_task
import "../tasks/concat_pon_vcf_task.wdl" as concat_pon_vcf_task

workflow PerbaseToPonVcfSplitChroms {
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

        String docker_image = "jiminpark/perbase-pon-vcf:1.4"
        Int threads = 2
        Int memSizeGB = 8
        Int diskSizeGB = 0

        Int index_memSizeGB = 4
        Int index_diskSizeGB = 0

        Int concat_memSizeGB = 4
        Int concat_diskSizeGB = 0
    }

    scatter (tsv in perbase_tsvs) {
        call index_perbase_tsv_task.index_perbase_tsv as index_tsv {
            input:
                tsv = tsv,
                docker_image = docker_image,
                memSizeGB = index_memSizeGB,
                diskSizeGB = if index_diskSizeGB > 0 then index_diskSizeGB else round(size(tsv, "GB")) + 10
        }
    }

    scatter (chrom in chroms) {
        call perbase_to_pon_vcf_task.perbase_to_pon_vcf as perbase_to_pon_vcf {
            input:
                perbase_tsvs = perbase_tsvs,
                perbase_tsv_indexes = index_tsv.tsv_index,
                samples = samples,
                reference_fasta = reference_fasta,
                reference_fasta_index = reference_fasta_index,
                pon_vcf_name = chrom + ".pon.vcf.gz",
                min_samples_covered = min_samples_covered,
                min_samples_alt = min_samples_alt,
                chroms = [chrom],
                header_chroms = chroms,
                docker_image = docker_image,
                threads = threads,
                memSizeGB = memSizeGB,
                diskSizeGB = if diskSizeGB > 0 then diskSizeGB else round(size(perbase_tsvs, "GB") + size(index_tsv.tsv_index, "GB")) + 20
        }
    }

    call concat_pon_vcf_task.concat_pon_vcf as concat_pon_vcf {
        input:
            vcfs = perbase_to_pon_vcf.pon_vcf,
            pon_vcf_name = pon_vcf_name,
            docker_image = docker_image,
            memSizeGB = concat_memSizeGB,
            diskSizeGB = if concat_diskSizeGB > 0 then concat_diskSizeGB else 2 * round(size(perbase_to_pon_vcf.pon_vcf, "GB")) + 10
    }

    output {
        File pon_vcf = concat_pon_vcf.pon_vcf
        File pon_vcf_index = concat_pon_vcf.pon_vcf_index
    }
}
