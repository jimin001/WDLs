version 1.0

import "../tasks/perbase_task.wdl" as perbase_task

workflow Perbase {
    input {
        File bam
        File bam_idx
        String sample

        Int min_base_quality_score = 20
        Int min_mapq = 10

        String docker_image = "jiminpark/perbase:1.4.0"
        Int threads = 8
        Int memSizeGB = 32
        Int diskSizeGB = 0
    }

    call perbase_task.perbase as perbase {
        input:
            bam = bam,
            bam_idx = bam_idx,
            sample = sample,
            min_base_quality_score = min_base_quality_score,
            min_mapq = min_mapq,
            docker_image = docker_image,
            threads = threads,
            memSizeGB = memSizeGB,
            diskSizeGB = if diskSizeGB > 0 then diskSizeGB else 3 * round(size(bam, "G")) + 50
    }

    output {
        File perbase_output = perbase.perbase_output
    }
}
