version 1.0

import "../tasks/check_local_noise.wdl" as noise

workflow CheckLocalNoise {
    input {
        File bam
        File bai
        File vcf
        File vcf_idx
        File ref
        File ref_fai
        File script

        String sample

        Boolean tag_vcf = false
        Int window = 10
        Int min_mapq = 20
    }

    call noise.check_local_noise {
        input:
            bam      = bam,
            bai      = bai,
            vcf      = vcf,
            vcf_idx  = vcf_idx,
            ref      = ref,
            ref_fai  = ref_fai,
            script   = script,
            sample   = sample,
            tag_vcf  = tag_vcf,
            window   = window,
            min_mapq = min_mapq
    }

    output {
        File        cleanliness_tsv = check_local_noise.cleanliness_tsv
        Array[File] tagged_vcf      = check_local_noise.tagged_vcf
        Array[File] tagged_vcf_idx  = check_local_noise.tagged_vcf_idx
    }
}
