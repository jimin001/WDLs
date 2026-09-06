version 1.0

# scores a pre-trained filter_nn model against a single sample via `filter_nn run`
# see: 08182026_deepsomatic_filter_nn_run_inference_5sample_16train_6test_pon_filtered_labeled_part2.sh

import "../tasks/DeepSomaticFilterNN_task.wdl" as filter_nn

workflow DeepSomaticFilterNN_run_inference {
    input {
        File model
        String model_type = "xgboost"
        Float threshold = 0.5

        String sample
        File vcf
        File vcf_idx
        File bam
        File bai

        # region to restrict scoring to, e.g. a held-out test bed
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

        String output_suffix = "output_labeled.wg.filtered"

        Int threads = 16
        Int memSizeGB = 300
    }

    call filter_nn.run_inference as run_inference {
        input:
            model = model,
            model_type = model_type,
            threshold = threshold,
            vcf = vcf,
            vcf_idx = vcf_idx,
            bam = bam,
            bai = bai,
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
            sample = sample,
            output_suffix = output_suffix,
            threads = threads,
            memSizeGB = memSizeGB
    }

    output {
        File filtered_vcf = run_inference.filtered_vcf
    }
}
