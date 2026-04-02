version 1.0

import "../tasks/merge_vcfs.wdl" as merge_vcfs
import "../tasks/DeepSomatic_tumoronly_task.wdl" as deepsomatic

# splits bam file to implement scatter-gather variant calling
workflow DeepSomatic_tumoronly_split_bams {
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


    output {
        File DeepSomaticVCF = merge_vcfs.merged_vcf
        File DeepSomaticVCFIndex = merge_vcfs.merged_vcf_idx
    }
}
