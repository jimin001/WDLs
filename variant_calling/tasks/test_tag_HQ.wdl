version 1.0

task tag_HQ {
    input {
        File bam
        File bai

        File germline_vcf
        File germline_vcf_idx
        File somatic_vcf
        File somatic_vcf_idx

        Int agreeing_gv_threshold = 5
        Int disagreeing_gv_threshold = 0
        Int unphased_threshold = 10
        String chromosome

        String vcf_base = basename(somatic_vcf, ".vcf.gz")
        
        String docker_image = "jiminpark/deepsomatic_postprocess:v3"
        Int threads = 4
        Int memSizeGB = 4
        Int diskSizeGB = round(size(bam, 'G')) * 2
    }

    command <<<
        # echo each line of the script to stdout so we can see what is happening
        # to turn off echo do 'set +o xtrace'
        set -o xtrace

        if [[ "~{chromosome}" == "" ]]
        then
                CHROMOSOME=""
        else
                CHROMOSOME="-c ~{chromosome}"
        fi
        
        mkdir -p "output_files"

        python3 /opt/scripts/tag_haplotag_exclusive_multiprocessing_v3.py \
        -bam ~{bam} \
        -v ~{somatic_vcf} \
        -g ~{germline_vcf} \
        -d "output_files" \
        -o "output_files/output_merged.vcf.gz" \
        -k ~{agreeing_gv_threshold} \
        -m ~{disagreeing_gv_threshold} \
        --max_workers ~{threads} \
        ${CHROMOSOME}

    >>>

    output {
        File output_vcf = "output_files/sorted_output_merged.vcf.gz"
        File output_vcf_idx = "output_files/sorted_output_merged.vcf.gz.tbi"

    }

    runtime {
        docker: docker_image
        cpu: threads
        memory: memSizeGB + " GB"
        disks: "local-disk " + diskSizeGB + " SSD"
    }
}