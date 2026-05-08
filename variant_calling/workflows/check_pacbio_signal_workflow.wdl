version 1.0

import "../tasks/check_pacbio_signal.wdl" as signal

workflow Check_pacbio_signal {
    input {
        File ont_vcf
        File ont_vcf_idx
        File pacbio_bam
        File pacbio_bai
        File script

        String sample

        String mode = "snv"
        Int min_mapq = 20
        Int min_baseq = 10
        Int min_alt_count = 1
        Float min_alt_freq = 0.0

        Int threads = 8
        Int memSizeGB = 16
        Int diskSizeGB = 0
    }

    call signal.check_pacbio_signal as check_pacbio_signal {
        input:
            ont_vcf = ont_vcf,
            ont_vcf_idx = ont_vcf_idx,
            pacbio_bam = pacbio_bam,
            pacbio_bai = pacbio_bai,
            script = script,
            sample = sample,
            mode = mode,
            min_mapq = min_mapq,
            min_baseq = min_baseq,
            min_alt_count = min_alt_count,
            min_alt_freq = min_alt_freq,
            threads = threads,
            memSizeGB = memSizeGB,
            diskSizeGB = if diskSizeGB > 0 then diskSizeGB else round(size(pacbio_bam, "G")) * 3
    }

    output {
        File signal_tsv = check_pacbio_signal.signal_tsv
    }
}
