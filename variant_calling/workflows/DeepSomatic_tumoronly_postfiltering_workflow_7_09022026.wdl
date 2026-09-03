version 1.0

import "../tasks/DeepSomatic_tumoronly_postfiltering_3.wdl" as filter

# this version of workflow:

workflow DeepSomatic_filtering {
    input {
        File somatic_VCF
        File somatic_IDX

        String sample

        File bam
        File bai

        Int agreeing_gv_threshold
        Int disagreeing_gv_threshold

        Boolean germline_is_harmphase = false
        File? germline_harmphase_VCF
        File? germline_harmphase_IDX
        File? germline_VCF
        File? germline_IDX

        File gnomad_VCF
        File gnomad_IDX
    }

    call filter.filter_pass_somatic as filter_pass_somatic {
        input: 
            somatic_VCF = somatic_VCF,
            somatic_IDX = somatic_IDX,
            sample = sample
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

    if (germline_is_harmphase) {
        call filter.extract_snps_from_harmphase as extract_snps_from_harmphase {
            input:
                harmphase_VCF = select_first([germline_harmphase_VCF]),
                sample = sample
        }
    }

    call filter.intersect_dv_and_gnomad as intersect_dv_and_gnomad {
        input:
            germline_VCF = select_first([extract_snps_from_harmphase.output_vcf, germline_VCF]),
            germline_IDX = select_first([extract_snps_from_harmphase.output_vcf_idx, germline_IDX]),
            gnomad_VCF = gnomad_VCF,
            gnomad_IDX = gnomad_IDX,
            sample = sample
    }

    call filter.tag_HQ as tag_HQ {
        input:
            bam = bam,
            bai = bai,
            germline_dv_intersect_gnomad = intersect_dv_and_gnomad.output_vcf,
            germline_dv_intersect_gnomad_idx = intersect_dv_and_gnomad.output_vcf_idx,
            somatic_VCF_input = subtract_segdup_regions.output_vcf,
            somatic_IDX_input = subtract_segdup_regions.output_vcf_idx
    }

    call filter.filter_hap_exclusive as filter_hap_exclusive {
        input:
            somatic_VCF_input = tag_HQ.output_vcf,
            somatic_IDX_input = tag_HQ.output_vcf_idx,
            sample = sample
    }

    call filter.filter_hap_vaf as filter_hap_vaf {
        input:
            somatic_VCF_input = filter_hap_exclusive.output_vcf,
            somatic_IDX_input = filter_hap_exclusive.output_vcf_idx,
            sample = sample
    }

    call filter.subtract_trhp_regions as subtract_trhp_regions {
        input:
            somatic_VCF_input = filter_hap_vaf.output_vcf,
            somatic_IDX_input = filter_hap_vaf.output_vcf_idx,
            sample = sample
    }

    output {
        File somatic_pass_only_vcf = filter_pass_somatic.somatic_pass_only_vcf
        File somatic_pass_only_vcf_idx = filter_pass_somatic.somatic_pass_only_vcf_idx

        File somatic_gnomad_tag = tag_gnomad_filter_all.gnomad_VCF
        File somatic_gnomad_tag_idx = tag_gnomad_filter_all.gnomad_IDX

        File somatic_gnomad_filter = tag_gnomad_filter_all.gnomad_filter_VCF
        File somatic_gnomad_filter_idx = tag_gnomad_filter_all.gnomad_filter_IDX

        File somatic_GQ20_DP10 = filter_GQ_and_DP.output_vcf
        File somatic_GQ20_DP10_idx = filter_GQ_and_DP.output_vcf_idx

        File somatic_segdup = subtract_segdup_regions.output_vcf
        File somatic_segdup_idx = subtract_segdup_regions.output_vcf_idx

        File soamtic_intersect_dv_and_gnomad = intersect_dv_and_gnomad.output_vcf
        File soamtic_intersect_dv_and_gnomad_idx = intersect_dv_and_gnomad.output_vcf_idx

        File somatic_tag_HQ_and_AH = tag_HQ.output_vcf
        File somatic_tag_HQ_and_AH_idx = tag_HQ.output_vcf_idx

        File somatic_hap_exclusive = filter_hap_exclusive.output_vcf
        File somatic_hap_exclusive_idx = filter_hap_exclusive.output_vcf_idx

        File somatic_hap_vaf = filter_hap_vaf.output_vcf
        File somatic_hap_vaf_idx = filter_hap_vaf.output_vcf_idx

        File somatic_trhp = subtract_trhp_regions.output_vcf
        File somatic_trhp_idx = subtract_trhp_regions.output_vcf_idx
    }
}