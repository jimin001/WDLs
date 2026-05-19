version 1.0

import "../tasks/DeepSomatic_tumoronly_postfiltering_2.wdl" as filter

workflow DeepSomatic_filtering {
    input {
        File somatic_VCF
        File somatic_IDX

        String sample

        File bam
        File bai

        File germline_dv_intersect_gnomad
        File germline_dv_intersect_gnomad_idx
    }

	call filter.tag_HQ as tag_HQ {
        input:
            bam = bam,
            bai = bai,
            germline_dv_intersect_gnomad = germline_dv_intersect_gnomad,
            germline_dv_intersect_gnomad_idx = germline_dv_intersect_gnomad_idx,
            somatic_VCF_input = somatic_VCF,
            somatic_IDX_input = somatic_IDX
    }

	output {
        File somatic_tag_HQ_and_AH = tag_HQ.output_vcf
        File somatic_tag_HQ_and_AH_idx = tag_HQ.output_vcf_idx
    }
}