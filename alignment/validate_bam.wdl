version 1.0

workflow validateBam {
    input {
        String sample_name
        File bam
        File? bam_index
    }

    call checkBamIntegrity {
        input:
            sample_name = sample_name,
            bam = bam,
            bam_index = bam_index
    }

    call averageCoverage {
        input:
            sample_name = sample_name,
            bam = bam
    }

    output {
        File validation_report = checkBamIntegrity.validation_report
        File flagstat = checkBamIntegrity.flagstat
        String total_reads = checkBamIntegrity.total_reads
        String mapped_reads = checkBamIntegrity.mapped_reads
        Float average_coverage = averageCoverage.average_coverage
        File coverage_report = averageCoverage.coverage_report
    }

}

task checkBamIntegrity {
    input {
        String docker_image = "jiminpark/aligners"

        String sample_name

        File bam
        File? bam_index

        Int threads = 8
        Int memSizeGB = 16
    }
    Int file_size = ceil(size(bam, "GB"))
    Int diskSizeGB = file_size + 20

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

        mkdir working
        cd working

        REPORT=~{sample_name}.validation_report.txt
        echo "sample: ~{sample_name}" > ${REPORT}
        echo "bam: ~{bam}" >> ${REPORT}

        # 1. quickcheck: valid header, and BGZF EOF marker present (catches files cut off at the end)
        samtools quickcheck -v ~{bam}
        echo "quickcheck: PASS" >> ${REPORT}

        # 2. decode every record; samtools errors out on corrupted/truncated blocks mid-file
        samtools flagstat -@ ~{threads} ~{bam} > ~{sample_name}.flagstat.txt
        echo "flagstat full read-through: PASS" >> ${REPORT}

        # printf %.0f so awk doesn't switch to scientific notation for billions of reads
        TOTAL=$(awk 'NR==1 {printf "%.0f\n", $1 + $3}' ~{sample_name}.flagstat.txt)
        MAPPED=$(grep -m1 " mapped (" ~{sample_name}.flagstat.txt | awk '{printf "%.0f\n", $1 + $3}')
        echo ${TOTAL} > total_reads.txt
        echo ${MAPPED} > mapped_reads.txt
        echo "total reads: ${TOTAL}" >> ${REPORT}
        echo "mapped reads: ${MAPPED}" >> ${REPORT}

        if [[ ${TOTAL} -eq 0 ]]
        then
            echo "ERROR: bam contains no reads" | tee -a ${REPORT} >&2
            exit 1
        fi

        # 3. bwa-mem2 output should be coordinate sorted
        if samtools view -H ~{bam} | grep -q "^@HD.*SO:coordinate"
        then
            echo "sort order: coordinate PASS" >> ${REPORT}
        else
            echo "ERROR: bam header is not SO:coordinate" | tee -a ${REPORT} >&2
            exit 1
        fi

        # record the program lines so we can confirm the aligner/version
        samtools view -H ~{bam} | grep "^@PG" >> ${REPORT} || echo "no @PG lines" >> ${REPORT}

        # 4. if an index was provided, make sure it matches the bam and is usable
        if [[ "~{bam_index}" != "" ]]
        then
            ln -s ~{bam} input.bam
            ln -s ~{bam_index} input.bam.bai
            samtools idxstats input.bam > ~{sample_name}.idxstats.txt
            echo "index (idxstats): PASS" >> ${REPORT}
        fi

        echo "OVERALL: PASS" >> ${REPORT}
    >>>

    output {
        File validation_report = "working/~{sample_name}.validation_report.txt"
        File flagstat = "working/~{sample_name}.flagstat.txt"
        # String, not Int: Cromwell's Int is 32-bit and WGS read counts can exceed 2^31
        String total_reads = read_string("working/total_reads.txt")
        String mapped_reads = read_string("working/mapped_reads.txt")
    }

    runtime {
        memory: memSizeGB + " GB"
        cpu: threads
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: docker_image
    }
}

task averageCoverage {
    input {
        String docker_image = "jiminpark/aligners"

        String sample_name

        File bam

        # false: average over covered positions only (plain `samtools depth`)
        # true: include zero-depth positions (`samtools depth -a`), i.e. genome-wide average
        Boolean include_zero_depth = false

        Int threads = 4
        Int memSizeGB = 8
    }
    Int file_size = ceil(size(bam, "GB"))
    Int diskSizeGB = file_size + 20

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

        mkdir working
        cd working

        if [[ "~{include_zero_depth}" == "true" ]]
        then
            DEPTH_ARGS="-a"
        else
            DEPTH_ARGS=""
        fi

        samtools depth -@ ~{threads} ${DEPTH_ARGS} ~{bam} \
        | awk '{sum+=$3} END { if (NR > 0) print sum/NR; else print 0 }' > average_coverage.txt

        echo "sample: ~{sample_name}" > ~{sample_name}.coverage.txt
        echo "include_zero_depth: ~{include_zero_depth}" >> ~{sample_name}.coverage.txt
        echo "Average = $(cat average_coverage.txt)" >> ~{sample_name}.coverage.txt
    >>>

    output {
        Float average_coverage = read_float("working/average_coverage.txt")
        File coverage_report = "working/~{sample_name}.coverage.txt"
    }

    runtime {
        memory: memSizeGB + " GB"
        cpu: threads
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: docker_image
    }
}
