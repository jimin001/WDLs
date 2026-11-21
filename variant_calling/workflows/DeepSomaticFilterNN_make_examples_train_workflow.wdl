version 1.0

# builds filter_nn training examples per sample (`filter_nn make-examples`), then trains one model
# on all of them (`filter_nn train`). Pair with DeepSomaticFilterNN_run_inference_workflow.wdl.
# see: slurm-jobs/SCRIPTS/card_analysis/10152026_deepsomatic_filter_nn_make_examples_train_11train11test_pon_filtered_vPoN2x2_NR5.sh

import "../tasks/DeepSomaticFilterNN_task.wdl" as filter_nn

struct FilterNNLabeledSample {
    String sample
    File bam
    File bai
    # labeled VCF: INFO/<label_field> = positive | negative
    File vcf
    File vcf_idx
}

workflow DeepSomaticFilterNN_make_examples_train {
    input {
        Array[FilterNNLabeledSample] training_samples

        # training regions shared by every sample, e.g. 11train_11test.bed
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

        # examples are named <sample>.<examples_suffix>.parquet
        String examples_suffix = "11train11test.train"
        String label_field = "type"
        Boolean methylation = true
        Int min_mapq = 20
        Int workers = 8

        # model is written as <model_prefix>.json (+ .columns.json, .importances.tsv)
        String model_prefix = "xgboost_somatic_filter"
        String model_type = "xgboost"

        # measured on the 11train_11test vPoN2x2_NR5 set: make-examples ~1 min / ~2 GB peak per
        # sample (reads at labeled sites are streamed from the BAM, not the whole BAM);
        # train 6 s / 0.24 GB peak
        Int make_examples_threads = 8
        Int make_examples_memSizeGB = 16
        Int train_threads = 4
        Int train_memSizeGB = 8
        Int preemptible = 3
    }

    scatter (training_sample in training_samples) {
        call filter_nn.make_examples as make_examples {
            input:
                bam = training_sample.bam,
                bai = training_sample.bai,
                vcf = training_sample.vcf,
                vcf_idx = training_sample.vcf_idx,
                region_bed = region_bed,
                reference = reference,
                reference_idx = reference_idx,
                reference_gzi = reference_gzi,
                gnomad_vcf = gnomad_vcf,
                gnomad_vcf_idx = gnomad_vcf_idx,
                cpg_island_bed = cpg_island_bed,
                cpg_island_bed_idx = cpg_island_bed_idx,
                cpg_shore_bed = cpg_shore_bed,
                cpg_shore_bed_idx = cpg_shore_bed_idx,
                phastcons_bed = phastcons_bed,
                phastcons_bed_idx = phastcons_bed_idx,
                segdup_bed = segdup_bed,
                segdup_bed_idx = segdup_bed_idx,
                sample = training_sample.sample,
                output_suffix = examples_suffix,
                label_field = label_field,
                methylation = methylation,
                min_mapq = min_mapq,
                workers = workers,
                threads = make_examples_threads,
                memSizeGB = make_examples_memSizeGB,
                preemptible = preemptible
        }
    }

    call filter_nn.train as train {
        input:
            train_parquet = make_examples.examples_parquet,
            output_prefix = model_prefix,
            model_type = model_type,
            threads = train_threads,
            memSizeGB = train_memSizeGB,
            preemptible = preemptible
    }

    output {
        Array[File] training_examples_parquet = make_examples.examples_parquet
        File model = train.model
        File model_columns = train.model_columns
        File model_importances = train.model_importances
    }
}
