version 1.0

workflow alignShortReadsOra {

    input {
        String sample_name

        File fastq1_ora
        File fastq2_ora

        # prebuilt bwa-mem2 index (see bwa_mem2_index.wdl); the fasta itself is not needed
        File ref_0123
        File ref_amb
        File ref_ann
        File ref_bwt_2bit_64
        File ref_pac

        String? additional_args

        # override the computed disk size for bwaMem2Ora
        Int? bwa_diskSizeGB
    }

    call bwaMem2Ora {
        input:
            sample_name = sample_name,
            fastq1_ora = fastq1_ora,
            fastq2_ora = fastq2_ora,
            ref_0123 = ref_0123,
            ref_amb = ref_amb,
            ref_ann = ref_ann,
            ref_bwt_2bit_64 = ref_bwt_2bit_64,
            ref_pac = ref_pac,
            additional_args = additional_args,
            diskSizeGB = bwa_diskSizeGB
    }

    output {
        File aligned_bam = bwaMem2Ora.aligned_bam
        File aligned_bam_idx = bwaMem2Ora.aligned_bam_idx
        File ora_summary = bwaMem2Ora.ora_summary
    }

}

task bwaMem2Ora {
    input {
        # built from Dockerfile.bwa_mem2_ora (jiminpark/aligners + orad, ORA_REF_PATH set)
        String docker_image = "jiminpark/bwa_mem2_ora:2.7.0"

        String sample_name

        File fastq1_ora
        File fastq2_ora

        File ref_0123
        File ref_amb
        File ref_ann
        File ref_bwt_2bit_64
        File ref_pac

        String? additional_args

        Int threads = 64
        Int orad_threads = 4
        Int sort_threads = 8
        Int memSizeGB = 96
        # long 64-core job is preempted often and restarts from scratch each time;
        # one cheap try, then on-demand
        Int preemptible = 1

        Int? diskSizeGB
    }
    # ORA is ~5x smaller than fastq.gz and the BAM is roughly fastq.gz-sized;
    # leave room for the ORA inputs + sorted BAM + samtools sort temp files + index
    Int ora_size = ceil(size(fastq1_ora, "GB") + size(fastq2_ora, "GB"))
    Int index_size = ceil(size([ref_0123, ref_amb, ref_ann, ref_bwt_2bit_64, ref_pac], "GB"))
    Int computedDiskSizeGB = 11 * ora_size + index_size + 50

    String out_bam = "~{sample_name}_Illumina.GRCh38.sorted.bam"

    command <<<
        # Set the exit code of a pipeline to that of the rightmost command
        # to exit with a non-zero status, or zero if all commands of the pipeline exit
        set -o pipefail
        # cause a bash script to exit immediately when a command fails
        set -e
        # cause the bash shell to treat unset variables as an error and exit immediately
        set -u
        # echo each line of the script to stdout so we can see what is happening
        # to turn off echo do 'set +o xtrace'
        set -o xtrace

        mkdir working ref
        cd working

        # check if length of "additional_args" is zero

        if [[ "~{additional_args}" == "" ]]
        then
            ADDITIONAL_ARGS=""
        else
            ADDITIONAL_ARGS="~{additional_args}"
        fi

        # index files may be localized to different directories; put them under one prefix
        ln -s ~{ref_0123} ../ref/ref.fa.0123
        ln -s ~{ref_amb} ../ref/ref.fa.amb
        ln -s ~{ref_ann} ../ref/ref.fa.ann
        ln -s ~{ref_bwt_2bit_64} ../ref/ref.fa.bwt.2bit.64
        ln -s ~{ref_pac} ../ref/ref.fa.pac

        # quick sanity check: ora file stats
        orad -i ~{fastq1_ora} > ora_summary.txt
        orad -i ~{fastq2_ora} >> ora_summary.txt

        # stream uncompressed fastq from orad straight into bwa-mem2 (no intermediate fastq.gz).
        # named pipes + explicit wait so an orad failure fails the task instead of
        # silently truncating the input, which <(...) process substitution would not catch
        mkfifo R1.fq R2.fq
        orad --raw --stdout -t ~{orad_threads} ~{fastq1_ora} > R1.fq &
        ORAD1_PID=$!
        orad --raw --stdout -t ~{orad_threads} ~{fastq2_ora} > R2.fq &
        ORAD2_PID=$!

        # align with bwa-mem2 and sort with samtools
        bwa-mem2 mem -M ${ADDITIONAL_ARGS} -t ~{threads} ../ref/ref.fa R1.fq R2.fq \
        | samtools sort -@ ~{sort_threads} -m 2G -o ~{out_bam} -

        # bwa-mem2 exits 0 on a truncated fastq, so the orad exit codes are the only check.
        # collect both before failing (if one orad dies the other can get SIGPIPE, 141)
        set +e
        wait ${ORAD1_PID}; ORAD1_RC=$?
        wait ${ORAD2_PID}; ORAD2_RC=$?
        set -e
        if [[ ${ORAD1_RC} -ne 0 || ${ORAD2_RC} -ne 0 ]]
        then
            echo "orad failed: R1 exit ${ORAD1_RC}, R2 exit ${ORAD2_RC}" >&2
            exit 1
        fi

        # index bam
        samtools index -@ ~{threads} ~{out_bam}
    >>>

    output {
        File aligned_bam = "working/~{out_bam}"
        File aligned_bam_idx = "working/~{out_bam}.bai"
        File ora_summary = "working/ora_summary.txt"
    }

    runtime {
        memory: memSizeGB + " GB"
        cpu: threads
        disks: "local-disk " + select_first([diskSizeGB, computedDiskSizeGB]) + " SSD"
        docker: docker_image
        preemptible: preemptible
    }
}
