version 1.0

import "../tasks/DeepSomatic_tumoronly_postfiltering_3.wdl" as filter

# this version of workflow: Illumina postfiltering
# PASS -> gnomAD filter -> GQ/DP -> segdup subtraction -> TR/HP subtraction
# (no haplotype-based filters, so no BAM or germline VCF needed)

workflow DeepSomatic_filtering {
    input {
        File somatic_VCF
        File somatic_IDX

        String sample
    }

    # downstream tasks name outputs off the input basename, so tagging here carries "Illumina" through every file
    String sample_tagged = sample + "_Illumina"

    call filter.filter_pass_somatic as filter_pass_somatic {
        input:
            somatic_VCF = somatic_VCF,
            somatic_IDX = somatic_IDX,
            sample = sample_tagged
    }

    call filter.tag_gnomad_filter_all as tag_gnomad_filter_all {
        input:
            somatic_VCF_input = filter_pass_somatic.somatic_pass_only_vcf,
            somatic_IDX_input = filter_pass_somatic.somatic_pass_only_vcf_idx,
            sample = sample
    }

    call filter.filter_GQ_and_DP as filter_GQ_and_DP {
        input:
            somatic_VCF_input = tag_gnomad_filter_all.gnomad_filter_VCF ,
            somatic_IDX_input = tag_gnomad_filter_all.gnomad_filter_IDX ,
            sample = sample
    }

    call filter.subtract_segdup_regions as subtract_segdup_regions {
        input:
            somatic_VCF_input = filter_GQ_and_DP.output_vcf,
            somatic_IDX_input = filter_GQ_and_DP.output_vcf_idx,
            sample = sample
    }

    call filter.subtract_trhp_regions as subtract_trhp_regions {
        input:
            somatic_VCF_input = subtract_segdup_regions.output_vcf,
            somatic_IDX_input = subtract_segdup_regions.output_vcf_idx,
            sample = sample
    }

    output {
        File somatic_pass_only_vcf = filter_pass_somatic.somatic_pass_only_vcf
        File somatic_pass_only_vcf_idx = filter_pass_somatic.somatic_pass_only_vcf_idx

        File somatic_gnomad_filter = tag_gnomad_filter_all.gnomad_filter_VCF
        File somatic_gnomad_filter_idx = tag_gnomad_filter_all.gnomad_filter_IDX

        File somatic_GQ20_DP10 = filter_GQ_and_DP.output_vcf
        File somatic_GQ20_DP10_idx = filter_GQ_and_DP.output_vcf_idx

        File somatic_segdup = subtract_segdup_regions.output_vcf
        File somatic_segdup_idx = subtract_segdup_regions.output_vcf_idx

        File somatic_trhp = subtract_trhp_regions.output_vcf
        File somatic_trhp_idx = subtract_trhp_regions.output_vcf_idx
    }
}
