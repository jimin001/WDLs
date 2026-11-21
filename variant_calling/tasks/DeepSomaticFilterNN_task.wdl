version 1.0

# google/deepvariant:filter_nn make-examples / train / run
# see: 07182026_deepsomatic_filter_nn_make_examples_RUSH_001_006_008_054_10mbp_regions.sh
#      07182026_deepsomatic_filter_nn_train_RUSH_001__006_008_054_10mbp_regions.sh
#      07182026_deepsomatic_filter_nn_run_inference_RUSH_034_labeled_wg_4sample_model.sh

task make_examples {
    input {
        File bam
        File bai
        File vcf
        File vcf_idx
        File region_bed

        File reference
        File reference_idx
        File reference_gzi

        File gnomad_vcf
        File gnomad_vcf_idx

        File cpg_island_bed
        File cpg_island_bed_idx
        File cpg_shore_bed
        File cpg_shore_bed_idx
        File phastcons_bed
        File phastcons_bed_idx
        File segdup_bed
        File segdup_bed_idx

        String sample
        String output_suffix = "train"
        String label_field = "type"
        Boolean methylation = true
        Int min_mapq = 20
        Int workers = 8

        String docker_image = "google/deepvariant:filter_nn"
        # filter_nn fetches reads chunk by chunk and evicts them, so it should not need much memory
        # (~2 GB peak on RUSH_008's labeled-site reads); 32 GB leaves room for the full BAM
        Int threads = 8
        Int memSizeGB = 32
        # the full BAM (~110 GB for ~30x ONT) is localized, plus gnomAD (~5 GB) and the reference
        Int diskSizeGB = ceil(size(bam, "GB") + size(gnomad_vcf, "GB") + size(reference, "GB")) + 50
        Int preemptible = 3
    }

    command <<<
        set -o pipefail
        set -e
        set -u
        set -o xtrace

        # The BAM is localized in full. Streaming only the labeled-site reads with
        # `samtools view -M -L gs://...` OOMs on Terra: htslib 1.16 (this image) leaks
        # ~30-45 MB/s while seeking over GCS, regardless of -@ (reproduced on MOSR_006, 523 sites).
        ln -s ~{bam} sample.bam
        ln -s ~{bai} sample.bam.bai
        ln -s ~{vcf} sample.vcf.gz
        ln -s ~{vcf_idx} sample.vcf.gz.tbi

        ln -s ~{reference} reference.fa.gz
        ln -s ~{reference_idx} reference.fa.gz.fai
        ln -s ~{reference_gzi} reference.fa.gz.gzi

        mkdir -p annotations
        ln -s ~{cpg_island_bed} annotations/cpg_islands.bed.gz
        ln -s ~{cpg_island_bed_idx} annotations/cpg_islands.bed.gz.tbi
        ln -s ~{cpg_shore_bed} annotations/cpg_shores.bed.gz
        ln -s ~{cpg_shore_bed_idx} annotations/cpg_shores.bed.gz.tbi
        ln -s ~{phastcons_bed} annotations/phastcons.bed.gz
        ln -s ~{phastcons_bed_idx} annotations/phastcons.bed.gz.tbi
        ln -s ~{segdup_bed} annotations/segdups.bed.gz
        ln -s ~{segdup_bed_idx} annotations/segdups.bed.gz.tbi
        ln -s ~{gnomad_vcf} annotations/gnomad.genomes.v4.1.sites.small.vcf.bgz
        ln -s ~{gnomad_vcf_idx} annotations/gnomad.genomes.v4.1.sites.small.vcf.bgz.tbi

        filter_nn make-examples \
            --bam sample.bam \
            --vcf sample.vcf.gz \
            --region ~{region_bed} \
            --ref reference.fa.gz \
            --output ~{sample}.~{output_suffix}.parquet \
            --label-field ~{label_field} \
            ~{true="--methylation" false="" methylation} \
            --min-mapq ~{min_mapq} \
            --workers ~{workers} \
            --bed-annotation cpg_island:annotations/cpg_islands.bed.gz:binary \
            --bed-annotation cpg_shore:annotations/cpg_shores.bed.gz:binary \
            --bed-annotation phastcons:annotations/phastcons.bed.gz:value \
            --bed-annotation segdup:annotations/segdups.bed.gz:binary \
            --gnomad annotations/gnomad.genomes.v4.1.sites.small.vcf.bgz
    >>>

    output {
        File examples_parquet = "~{sample}.~{output_suffix}.parquet"
    }

    runtime {
        preemptible: preemptible
        docker: docker_image
        cpu: threads
        memory: memSizeGB + " GB"
        disks: "local-disk " + diskSizeGB + " HDD"
    }
}

task train {
    input {
        Array[File] train_parquet

        String output_prefix
        String model_type = "xgboost"

        String docker_image = "google/deepvariant:filter_nn"
        # xgboost on 5 samples x ~500-750 rows: 6 s, 0.24 GB peak on 4 CPUs
        Int threads = 4
        Int memSizeGB = 8
        Int diskSizeGB = ceil(size(train_parquet, "GB")) + 20
        Int preemptible = 3
    }

    command <<<
        set -o pipefail
        set -e
        set -u
        set -o xtrace

        filter_nn train \
            --input ~{sep=" --input " train_parquet} \
            --output ~{output_prefix}.json \
            --model-type ~{model_type}
    >>>

    output {
        File model = "~{output_prefix}.json"
        File model_columns = "~{output_prefix}.json.columns.json"
        File model_importances = "~{output_prefix}.json.importances.tsv"
    }

    runtime {
        preemptible: preemptible
        docker: docker_image
        cpu: threads
        memory: memSizeGB + " GB"
        disks: "local-disk " + diskSizeGB + " HDD"
    }
}

task run_inference {
    input {
        File model
        String model_type = "xgboost"
        Float threshold = 0.5

        File vcf
        File vcf_idx
        File bam
        File bai
        File? region_bed

        File reference
        File reference_idx
        File reference_gzi

        File gnomad_vcf
        File gnomad_vcf_idx

        File cpg_island_bed
        File cpg_island_bed_idx
        File cpg_shore_bed
        File cpg_shore_bed_idx
        File phastcons_bed
        File phastcons_bed_idx
        File segdup_bed
        File segdup_bed_idx

        String sample
        String output_suffix = "output_labeled.wg.filtered"

        String docker_image = "google/deepvariant:filter_nn"
        Int threads = 16
        Int memSizeGB = 300
        Int diskSizeGB = 5 * round(size(bam, "G")) + 50
    }

    command <<<
        set -o pipefail
        set -e
        set -u
        set -o xtrace

        ln -s ~{bam} sample.bam
        ln -s ~{bai} sample.bam.bai
        ln -s ~{vcf} sample.vcf.gz
        ln -s ~{vcf_idx} sample.vcf.gz.tbi

        ln -s ~{reference} reference.fa.gz
        ln -s ~{reference_idx} reference.fa.gz.fai
        ln -s ~{reference_gzi} reference.fa.gz.gzi

        mkdir -p annotations
        ln -s ~{cpg_island_bed} annotations/cpg_islands.bed.gz
        ln -s ~{cpg_island_bed_idx} annotations/cpg_islands.bed.gz.tbi
        ln -s ~{cpg_shore_bed} annotations/cpg_shores.bed.gz
        ln -s ~{cpg_shore_bed_idx} annotations/cpg_shores.bed.gz.tbi
        ln -s ~{phastcons_bed} annotations/phastcons.bed.gz
        ln -s ~{phastcons_bed_idx} annotations/phastcons.bed.gz.tbi
        ln -s ~{segdup_bed} annotations/segdups.bed.gz
        ln -s ~{segdup_bed_idx} annotations/segdups.bed.gz.tbi
        ln -s ~{gnomad_vcf} annotations/gnomad.genomes.v4.1.sites.small.vcf.bgz
        ln -s ~{gnomad_vcf_idx} annotations/gnomad.genomes.v4.1.sites.small.vcf.bgz.tbi

        filter_nn run \
            --model ~{model} \
            --model-type ~{model_type} \
            --threshold ~{threshold} \
            ~{"--region " + region_bed} \
            --vcf sample.vcf.gz \
            --bam sample.bam \
            --ref reference.fa.gz \
            --bed-annotation cpg_island:annotations/cpg_islands.bed.gz:binary \
            --bed-annotation cpg_shore:annotations/cpg_shores.bed.gz:binary \
            --bed-annotation phastcons:annotations/phastcons.bed.gz:value \
            --bed-annotation segdup:annotations/segdups.bed.gz:binary \
            --gnomad annotations/gnomad.genomes.v4.1.sites.small.vcf.bgz \
            --output ~{sample}_~{output_suffix}.vcf.gz
    >>>

    output {
        File filtered_vcf = "~{sample}_~{output_suffix}.vcf.gz"
    }

    runtime {
        docker: docker_image
        cpu: threads
        memory: memSizeGB + " GB"
        disks: "local-disk " + diskSizeGB + " SSD"
    }
}
