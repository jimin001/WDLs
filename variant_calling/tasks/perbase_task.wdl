version 1.0

task perbase {
    input {
        File bam
        File bam_idx
        String sample

        File? bed_file
        Int min_base_quality_score = 10
        Int min_mapq = 10
        # 3848 = secondary (256) + QC-fail (512) + duplicate (1024) + supplementary (2048)
        Int exclude_flags = 3848
        Int compression_level = 6

        String docker_image = "jiminpark/perbase:1.4.0-bgzip"
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

        # perbase writes 20 columns; keep only:
        # REF POS DEPTH A C G T INS DEL REF_SKIP FAIL NEAR_MAX_DEPTH
        # (drops N, the IUPAC R/Y/S/W/K/M columns, and COUNT_OF_MATE_RESOLUTIONS)
        perbase base-depth ~{sample}.bam \
            -t ~{threads} \
            -F ~{exclude_flags} \
            --min-base-quality-score ~{min_base_quality_score} \
            --min-mapq ~{min_mapq} \
            ~{"-b " + bed_file} \
            | cut -f 1-7,15-18,20 \
            | bgzip -@ ~{threads} -l ~{compression_level} \
            > ~{sample}_perbase_output.tsv.gz
    >>>

    output {
        File perbase_output = "~{sample}_perbase_output.tsv.gz"
    }

    runtime {
        preemptible: 2
        docker: docker_image
        cpu: threads
        memory: memSizeGB + " GB"
        disks: "local-disk " + diskSizeGB + " SSD"
    }
}
