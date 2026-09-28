version 1.0

# Tabix-indexes a bgzipped perbase TSV (the Perbase WDL's
# *_perbase_output.tsv.gz) so perbase_to_pon_vcf can fetch one chromosome
# directly instead of streaming the whole file. Building the index also
# checks sortedness: tabix refuses positions out of order within a contig or
# a contig split into non-contiguous blocks, so this replaces the
# sort_perbase_tsv check for the per-chromosome workflow.
#
# Only the .tbi is output (not a copy of the TSV); pass it alongside the
# original TSV downstream.
task index_perbase_tsv {
    input {
        File tsv

        String docker_image = "jiminpark/perbase-pon-vcf:1.4"
        Int threads = 1
        Int memSizeGB = 4
        Int diskSizeGB = round(size(tsv, "GB")) + 10
    }

    command <<<
        set -o pipefail
        set -e
        set -u
        set -o xtrace

        # tabix needs BGZF: gzip magic + FEXTRA flag (1f 8b 08 04) and the
        # "BC" extra subfield id at bytes 13-14. Plain or plain-gzip TSVs
        # (older perbase output) must be bgzipped first.
        magic=$(head -c 14 ~{tsv} | od -An -tx1 | tr -d ' \n')
        if [ "${magic:0:8}" != "1f8b0804" ] || [ "${magic:24:4}" != "4243" ]; then
            echo "ERROR: ~{tsv} is not bgzip-compressed; recompress it with bgzip before indexing" >&2
            exit 1
        fi

        # tabix writes <file>.tbi next to the file it indexes, so index a
        # link in the working directory rather than the localized input.
        ln -s ~{tsv} perbase.tsv.gz

        python3 - <<'EOF'
import gzip
import pysam

with gzip.open("perbase.tsv.gz", "rt") as fh:
    header = fh.readline().rstrip("\n").split("\t")
if header[:2] != ["REF", "POS"]:
    raise SystemExit(f"ERROR: expected perbase's REF and POS as the first two columns; got {header[:2]}")

# REF (contig) and POS columns, 1-based, skipping the header line
pysam.tabix_index("perbase.tsv.gz", seq_col=0, start_col=1, end_col=1,
                  line_skip=1, zerobased=False, force=True)
EOF
    >>>

    output {
        File tsv_index = "perbase.tsv.gz.tbi"
    }

    runtime {
        preemptible: 2
        docker: docker_image
        cpu: threads
        memory: memSizeGB + " GB"
        disks: "local-disk " + diskSizeGB + " SSD"
    }
}
