version 1.0

import "../tasks/merge_vcfs.wdl" as merge_vcfs
import "../tasks/DeepSomatic_tumoronly_postfiltering.wdl" as filter
import "../tasks/DeepVariant_task.wdl" as deepvariant


workflow DeepVariant_split_bams {
    input {
        File bam
        File? bai

        Array[String] regions

        String output_prefix
        String sample
    }

    scatter (region in regions) {
        call deepvariant.split_bams as split_bams {
            input:
                    bam = bam,
                    bai = bai,
                    region = region,
                    sample = sample
        }

        call deepvariant.deepVariant as deepVariant {
            input:
                    tumor_bam = split_bams.region_bam,
                    tumor_bam_idx = split_bams.region_bai,
                    output_prefix = output_prefix
        }
    }

    call merge_vcfs.merge_vcfs as merge_vcfs {
        input:
                vcf_files = deepVariant.output_vcf,
                output_prefix = output_prefix
    }

    call filter.filter_pass_germline as filter_pass_germline {
        input: 
            germline_VCF = merge_vcfs.merged_vcf,
            sample = sample
    }


    output {
        File DeepVariantVCF = merge_vcfs.merged_vcf
        File DeepVariantVCFIndex = merge_vcfs.merged_vcf_idx

        File germline_pass_only_vcf = filter_pass_germline.germline_pass_only_vcf
        File germline_pass_only_vcf_idx = filter_pass_germline.germline_pass_only_vcf_idx
    }
}
