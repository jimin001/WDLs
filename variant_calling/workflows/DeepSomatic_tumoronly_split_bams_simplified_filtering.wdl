version 1.0

import "../tasks/merge_vcfs.wdl" as merge_vcfs
import "../tasks/DeepSomatic_tumoronly_postfiltering_2.wdl" as filter
import "../tasks/DeepSomatic_tumoronly_task.wdl" as deepsomatic

# 1/30/2026 simplified filtering steps
workflow DeepSomatic_tumoronly_split_bams_simplified_filtering {
    input {
        File bam
        File? bai

        Array[String] regions

        #File germline_vcf

        String output_prefix
        String sample
    }

    scatter (region in regions) {
        call deepsomatic.split_bams as split_bams {
            input:
                    bam = bam,
                    bai = bai,
                    region = region,
                    sample = sample
        }

        call deepsomatic.deepSomatic as deepSomatic {
            input:
                    tumor_bam = split_bams.region_bam,
                    tumor_bam_idx = split_bams.region_bai,
                    output_prefix = output_prefix,
                    sample_name_tumor = sample
        }
    }

    call merge_vcfs.merge_vcfs as merge_vcfs {
        input:
                vcf_files = deepSomatic.output_vcf,
                output_prefix = output_prefix
    }

    call filter.filter_pass_somatic as filter_pass_somatic {
        input: 
            somatic_VCF = merge_vcfs.merged_vcf,
            somatic_IDX = merge_vcfs.merged_vcf_idx,
            sample = sample
    }

    call filter.tag_gnomad as tag_gnomad {
        input:
            somatic_VCF_input = filter_pass_somatic.somatic_pass_only_vcf,
            somatic_IDX_input = filter_pass_somatic.somatic_pass_only_vcf_idx,
            sample = sample
    }

    call filter.filter_GQ_and_DP as filter_GQ_and_DP {
        input:
            somatic_VCF_input = tag_gnomad.gnomad_filter_VCF ,
            somatic_IDX_input = tag_gnomad.gnomad_filter_IDX ,
            sample = sample
    }

    call filter.subtract_segdup_regions as subtract_segdup_regions {
        input:
            somatic_VCF_input = filter_GQ_and_DP.output_vcf,
            somatic_IDX_input = filter_GQ_and_DP.output_vcf_idx,
            sample = sample
    }

    call filter.filter_high_vaf as filter_high_vaf {
        input:
            somatic_VCF_input = subtract_segdup_regions.output_vcf,
            somatic_IDX_input = subtract_segdup_regions.output_vcf_idx,
            sample = sample
    }

    output {
        File DeepSomaticVCF = merge_vcfs.merged_vcf
        File DeepSomaticVCFIndex = merge_vcfs.merged_vcf_idx

        File somatic_pass_only_vcf = filter_pass_somatic.somatic_pass_only_vcf
        File somatic_pass_only_vcf_idx = filter_pass_somatic.somatic_pass_only_vcf_idx

        File somatic_gnomad_tag = tag_gnomad.gnomad_VCF
        File somatic_gnomad_tag_idx = tag_gnomad.gnomad_IDX

        File somatic_gnomad_filter = tag_gnomad.gnomad_filter_VCF
        File somatic_gnomad_filter_idx = tag_gnomad.gnomad_filter_IDX

        File somatic_GQ20_DP10 = filter_GQ_and_DP.output_vcf
        File somatic_GQ20_DP10_idx = filter_GQ_and_DP.output_vcf_idx

        File somatic_segdup = subtract_segdup_regions.output_vcf
        File somatic_segdup_idx = subtract_segdup_regions.output_vcf_idx

        File somatic_filter_high_vaf = filter_high_vaf.output_vcf
        File somatic_filter_high_vaf_idx = filter_high_vaf.output_vcf_idx
    }
}
