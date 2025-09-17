version 1.0

import "../tasks/DeepSomatic_tumoronly_postfiltering.wdl" as filter

workflow DeepSomatic_filtering {
    input {
        File germline_VCF
        File germline_IDX

        File somatic_VCF
        File somatic_IDX

        String sample

        File bam
        File bai
    }

    call filter.filter_pass as filter_pass {
        input: 
            germline_VCF = germline_VCF,
            germline_IDX = germline_IDX,
            somatic_VCF = somatic_VCF,
            somatic_IDX = somatic_IDX,
            sample = sample
    }

    call filter.filter_out_germline_variants as filter_out_germline_variants {
        input:
            germline_VCF_pass_only = filter_pass.germline_pass_only_vcf,
            germline_IDX_pass_only = filter_pass.germline_pass_only_vcf_idx,
            somatic_VCF_pass_only = filter_pass.somatic_pass_only_vcf,
            somatic_IDX_pass_only = filter_pass.somatic_pass_only_vcf_idx,
            sample = sample
    }

    call filter.filter_GQ_and_DP as filter_GQ_and_DP {
        input:
            deepsomatic_only_VCF = filter_out_germline_variants.deepsomatic_only_vcf,
            deepsomatic_only_IDX = filter_out_germline_variants.deepsomatic_only_vcf_idx,
            sample = sample
    }

    call filter.subtract_segdup_regions as subtract_segdup_regions {
        input:
            deepsomatic_only_GQ20_DP10_VCF = filter_GQ_and_DP.deepsomatic_only_GQ20_DP10_vcf,
            deepsomatic_only_GQ20_DP10_IDX = filter_GQ_and_DP.deepsomatic_only_GQ20_DP10_vcf_idx,
            sample = sample
    }

    call filter.tag_haplotype as tag_haplotype {
        input:
            deepsomatic_only_GQ20_DP10_segdup_VCF = subtract_segdup_regions.deepsomatic_only_GQ20_DP10_segdup_vcf,
            deepsomatic_only_GQ20_DP10_segdup_IDX = subtract_segdup_regions.deepsomatic_only_GQ20_DP10_segdup_idx,
            BAM = bam,
            BAI = bai,
            sample = sample
    }

    call filter.tag_gnomad as tag_gnomad {
        input:
            deepsomatic_only_GQ20_DP10_segdup_tagAH_VCF = tag_haplotype.deepsomatic_only_GQ20_DP10_segdup_tagAH_vcf,
            deepsomatic_only_GQ20_DP10_segdup_tagAH_IDX = tag_haplotype.deepsomatic_only_GQ20_DP10_segdup_tagAH_idx,
            sample = sample
    }

    output {
        File germline_pass_only_vcf = filter_pass.germline_pass_only_vcf
        File germline_pass_only_vcf_idx = filter_pass.germline_pass_only_vcf_idx
        File somatic_pass_only_vcf = filter_pass.somatic_pass_only_vcf
        File somatic_pass_only_vcf_idx = filter_pass.somatic_pass_only_vcf_idx

        File deepsomatic_only_vcf = filter_out_germline_variants.deepsomatic_only_vcf
        File deepsomatic_only_vcf_idx = filter_out_germline_variants.deepsomatic_only_vcf_idx

        File deepsomatic_only_GQ20_DP10_vcf = filter_GQ_and_DP.deepsomatic_only_GQ20_DP10_vcf
        File deepsomatic_only_GQ20_DP10_vcf_idx = filter_GQ_and_DP.deepsomatic_only_GQ20_DP10_vcf_idx

        File deepsomatic_only_GQ20_DP10_segdup_vcf = subtract_segdup_regions.deepsomatic_only_GQ20_DP10_segdup_vcf
        File deepsomatic_only_GQ20_DP10_segdup_idx = subtract_segdup_regions.deepsomatic_only_GQ20_DP10_segdup_idx

        File deepsomatic_only_GQ20_DP10_segdup_tagAH_vcf = tag_haplotype.deepsomatic_only_GQ20_DP10_segdup_tagAH_vcf
        File deepsomatic_only_GQ20_DP10_segdup_tagAH_idx = tag_haplotype.deepsomatic_only_GQ20_DP10_segdup_tagAH_idx

        File gnomad_VCF = tag_gnomad.gnomad_VCF
        File gnomad_IDX = tag_gnomad.gnomad_IDX
        File gnomad_filter_VCF = tag_gnomad.gnomad_filter_VCF
        File gnomad_filter_IDX = tag_gnomad.gnomad_filter_IDX
    }
}