version 1.0

# Joins per-chromosome PoN VCF shards from perbase_to_pon_vcf into one
# bgzipped+tabix-indexed VCF. The shards must all have been written with the
# same --header-chroms (so their headers are identical) and be given in
# contig order; that lets bcftools concat --naive join their BGZF blocks
# as-is without decompressing or re-parsing any records.
task concat_pon_vcf {
    input {
        Array[File] vcfs
        String pon_vcf_name = "pon.vcf.gz"

        String docker_image = "jiminpark/perbase-pon-vcf:1.4"
        Int threads = 1
        Int memSizeGB = 4
        Int diskSizeGB = 2 * round(size(vcfs, "GB")) + 10
    }

    command <<<
        set -o pipefail
        set -e
        set -u
        set -o xtrace

        python3 - ~{pon_vcf_name} ~{sep=' ' vcfs} <<'EOF'
import sys
import pysam
import pysam.bcftools

out, shards = sys.argv[1], sys.argv[2:]
pysam.bcftools.concat("--naive", "-o", out, *shards, catch_stdout=False)
pysam.tabix_index(out, preset="vcf", force=True)
EOF
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
