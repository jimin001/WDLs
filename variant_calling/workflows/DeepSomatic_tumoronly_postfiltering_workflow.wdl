version 1.0

import ../tasks/DeepSomatic_tumoronly_postfiltering.wdl as filter

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

    call filter.filter_pass {
        input: 
            germline_VCF = germline_VCF,
            germline_IDX = germline_IDX,
            somatic_VCF = somatic_VCF,
            somatic_IDX = somatic_IDX,
            sample = sample
    }

    call filter.filter_out_germline_variants {
        input:
            germline_VCF_pass_only = filter.filter_pass.germline_pass_only_vcf,
            germline_IDX_pass_only = filter.filter_pass.germline_pass_only_vcf_idx,
            somatic_VCF_pass_only = filter.filter_pass.somatic_pass_only_vcf,
            somatic_IDX_pass_only = filter.filter_pass.somatic_pass_only_vcf_idx,
            sample = sample

    }

    call filter.filter_GQ_and_DP {
        input:
            deepsomatic_only_VCF = filter.filter_out_germline_variants.deepsomatic_only_vcf,
            deepsomatic_only_IDX = filter.filter_out_germline_variants.deepsomatic_only_idx,
            sample = sample
    }

    call filter.subtract_segdup_regions {
        input:
            deepsomatic_only_GQ20_DP10_VCF = filter.filter_GQ_and_DP.omatic_only_GQ20_DP10_vcf,
            deepsomatic_only_GQ20_DP10_IDX = filter.filter_GQ_and_DP.omatic_only_GQ20_DP10_vcf_idx,
            sample = sample
    }

    call filter.tag_haplotype {
        input:
            deepsomatic_only_GQ20_DP10_segdup_VCF = filter.subtract_segdup_regions.deepsomatic_only_GQ20_DP10_segdup_vcf,
            deepsomatic_only_GQ20_DP10_segdup_IDX = filter.subtract_segdup_regions.deepsomatic_only_GQ20_DP10_segdup_idx,
            BAM = bam,
            BAI = bai,
            sample = sample
    }

    call filter.tag_gnomad {
        input:
            deepsomatic_only_GQ20_DP10_segdup_tagAH_VCF = filter.tag_haplotype.deepsomatic_only_GQ20_DP10_segdup_tagAH_vcf,
            deepsomatic_only_GQ20_DP10_segdup_tagAH_IDX = filter.tag_haplotype.deepsomatic_only_GQ20_DP10_segdup_tagAH_idx,
            sample = sample
    }
}