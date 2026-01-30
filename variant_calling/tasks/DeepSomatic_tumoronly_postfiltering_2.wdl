version 1.0

################### tasks #####################

task filter_pass {
    input {
        File germline_VCF
        File somatic_VCF
        File somatic_IDX

        String sample

        String docker_image = "jiminpark/deepsomatic_postprocess"
        Int threads = 1
        Int memSizeGB = 4
        Int diskSizeGB = 64
    }

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

        bcftools index -t ~{germline_VCF}

        bcftools filter -i 'FILTER="PASS"' ~{germline_VCF} | bgzip > ~{sample}_deepvariant_pass_only.vcf.gz
        bcftools filter -i 'FILTER="PASS"' ~{somatic_VCF} | bgzip > ~{sample}_deepsomatic_pass_only.vcf.gz

        bcftools index -t ~{sample}_deepvariant_pass_only.vcf.gz
        bcftools index -t ~{sample}_deepsomatic_pass_only.vcf.gz
    >>>

    output {
        File germline_pass_only_vcf = "~{sample}_deepvariant_pass_only.vcf.gz"
        File germline_pass_only_vcf_idx = "~{sample}_deepvariant_pass_only.vcf.gz.tbi"
        File somatic_pass_only_vcf = "~{sample}_deepsomatic_pass_only.vcf.gz"
        File somatic_pass_only_vcf_idx = "~{sample}_deepsomatic_pass_only.vcf.gz.tbi"
        
    }

    runtime {
        docker: docker_image
        cpu: threads
        memory: memSizeGB + " GB"
        disks: "local-disk " + diskSizeGB + " SSD"
    }
}

task filter_pass_somatic {
    input {
        File somatic_VCF
        File somatic_IDX

        String sample

        String docker_image = "jiminpark/deepsomatic_postprocess"
        Int threads = 1
        Int memSizeGB = 4
        Int diskSizeGB = 64
    }

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

        bcftools filter -i 'FILTER="PASS"' ~{somatic_VCF} | bgzip > ~{sample}_deepsomatic_pass_only.vcf.gz

        bcftools index -t ~{sample}_deepsomatic_pass_only.vcf.gz
    >>>

    output {
        File somatic_pass_only_vcf = "~{sample}_deepsomatic_pass_only.vcf.gz"
        File somatic_pass_only_vcf_idx = "~{sample}_deepsomatic_pass_only.vcf.gz.tbi"
        
    }

    runtime {
        docker: docker_image
        cpu: threads
        memory: memSizeGB + " GB"
        disks: "local-disk " + diskSizeGB + " SSD"
    }
}

# updated task so it fits better as first step in workflow
task tag_gnomad {
    input {
        File somatic_VCF_input
        File somatic_IDX_input

        String sample
        String vcf_base = basename(somatic_VCF_input, ".vcf.gz")

        String docker_image = "jiminpark/snpeff"
        Int threads = 1
        Int memSizeGB = 16
        Int diskSizeGB = 64
    }

    command <<<
        # echo each line of the script to stdout so we can see what is happening
        # to turn off echo do 'set +o xtrace'
        set -o xtrace


        java -jar /home/apps/snpEff/SnpSift.jar Annotate -noId /opt/scripts/gnomad.genomes.v4.1.sites.small.vcf.bgz \
        ~{somatic_VCF_input} | bgzip > ~{vcf_base}_gnomad.vcf.gz

        bcftools index -t ~{vcf_base}_gnomad.vcf.gz

        #bcftools view -i 'INFO/AF<0.001 | INFO/AF="."' ~{vcf_base}_gnomad.vcf.gz | bgzip  > ~{vcf_base}_gnomadrare.vcf.gz
        bcftools view -i 'INFO/AF<0.0005 | INFO/AF="."' ~{vcf_base}_gnomad.vcf.gz | bgzip  > ~{vcf_base}_gnomadrare.vcf.gz

        bcftools index -t ~{vcf_base}_gnomadrare.vcf.gz

    >>>

    output {
        File gnomad_VCF = "~{vcf_base}_gnomad.vcf.gz"
        File gnomad_IDX = "~{vcf_base}_gnomad.vcf.gz.tbi"
        File gnomad_filter_VCF = "~{vcf_base}_gnomadrare.vcf.gz"
        File gnomad_filter_IDX = "~{vcf_base}_gnomadrare.vcf.gz.tbi"
    }

    runtime {
        docker: docker_image
        cpu: threads
        memory: memSizeGB + " GB"
        disks: "local-disk " + diskSizeGB + " SSD"
    }
}
task filter_pass_germline {
    input {
        File germline_VCF

        String sample

        String docker_image = "jiminpark/deepsomatic_postprocess"
        Int threads = 1
        Int memSizeGB = 4
        Int diskSizeGB = 64
    }

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

        bcftools index -t ~{germline_VCF}

        bcftools filter -i 'FILTER="PASS"' ~{germline_VCF} | bgzip > ~{sample}_deepvariant_pass_only.vcf.gz

        bcftools index -t ~{sample}_deepvariant_pass_only.vcf.gz
    >>>

    output {
        File germline_pass_only_vcf = "~{sample}_deepvariant_pass_only.vcf.gz"
        File germline_pass_only_vcf_idx = "~{sample}_deepvariant_pass_only.vcf.gz.tbi"        
    }

    runtime {
        docker: docker_image
        cpu: threads
        memory: memSizeGB + " GB"
        disks: "local-disk " + diskSizeGB + " SSD"
    }
}


task filter_GQ_and_DP {
    input {
        File somatic_VCF_input
        File somatic_IDX_input

        String sample
        String vcf_base = basename(somatic_VCF_input, ".vcf.gz")

        String docker_image = "jiminpark/deepsomatic_postprocess"
        Int threads = 1
        Int memSizeGB = 4
        Int diskSizeGB = 64
    }

    command <<<
        # echo each line of the script to stdout so we can see what is happening
        # to turn off echo do 'set +o xtrace'
        set -o xtrace

        bcftools view -e 'FMT/GQ<20' ~{somatic_VCF_input} | bcftools view -e 'FMT/DP<10' | bgzip > ~{vcf_base}_GQ20_DP10.vcf.gz
        bcftools index -t ~{vcf_base}_GQ20_DP10.vcf.gz
    >>>

    output {
        File output_vcf = "~{vcf_base}_GQ20_DP10.vcf.gz"
        File output_vcf_idx = "~{vcf_base}_GQ20_DP10.vcf.gz.tbi"
    }

    runtime {
        docker: docker_image
        cpu: threads
        memory: memSizeGB + " GB"
        disks: "local-disk " + diskSizeGB + " SSD"
    }
}

task subtract_segdup_regions {
    input {
        File somatic_VCF_input
        File somatic_IDX_input

        String sample
        String vcf_base = basename(somatic_VCF_input, ".vcf.gz")

        String docker_image = "jiminpark/deepsomatic_postprocess"
        Int threads = 1
        Int memSizeGB = 4
        Int diskSizeGB = 64
    }

    command <<<
        # echo each line of the script to stdout so we can see what is happening
        # to turn off echo do 'set +o xtrace'
        set -o xtrace

        bedtools subtract -header -a ~{somatic_VCF_input} -b /opt/scripts/GRCh38_segdups.bed | bgzip > ~{vcf_base}_segdup.vcf.gz
        bcftools index -t ~{vcf_base}_segdup.vcf.gz
    >>>

    output {
        File output_vcf = "~{vcf_base}_segdup.vcf.gz"
        File output_vcf_idx = "~{vcf_base}_segdup.vcf.gz.tbi"
    }

    runtime {
        docker: docker_image
        cpu: threads
        memory: memSizeGB + " GB"
        disks: "local-disk " + diskSizeGB + " SSD"
    }
}

task filter_high_vaf {
    input {
        File somatic_VCF_input
        File somatic_IDX_input

        String sample
        String vcf_base = basename(somatic_VCF_input, ".vcf.gz")

        String docker_image = "jiminpark/deepsomatic_postprocess"
        Int threads = 1
        Int memSizeGB = 4
        Int diskSizeGB = 64
    }

    command <<<
        set -o xtrace
        bcftools view -i 'FMT/VAF<=0.8' ~{somatic_VCF_input} | bgzip > ~{vcf_base}_VAF0.8.vcf.gz
        bcftools index -t  ~{vcf_base}_VAF0.8.vcf.gz
    >>>

    output {
        File output_vcf = "~{vcf_base}_VAF0.8.vcf.gz"
        File output_vcf_idx = "~{vcf_base}_VAF0.8.vcf.gz.tbi"
    }
}



task tag_haplotype {
    input {
        File somatic_VCF_input
        File somatic_IDX_input
        File BAM
        File? BAI
        
        String sample
        String vcf_base = basename(somatic_VCF_input, ".vcf.gz")

        String docker_image = "jiminpark/deepsomatic_postprocess"
        Int threads = 4
        Int memSizeGB = 4
        Int diskSizeGB = round(size(BAM, 'G')) + 100
    }

    command <<<
        if [[ "~{BAI}" == "" ]]
        then
                samtools index -@ ~{threads} ~{BAM}
        else
                continue
        fi

        # echo each line of the script to stdout so we can see what is happening
        # to turn off echo do 'set +o xtrace'
        set -o xtrace

        python3 /opt/scripts/haplotype_fraction.py -v ~{somatic_VCF_input} -bam ~{BAM} -o ~{vcf_base}_tagAH.vcf -u 10
        bgzip ~{vcf_base}_tagAH.vcf
        bcftools index -t ~{vcf_base}_tagAH.vcf.gz
    >>>

    output {
        File output_vcf = "~{vcf_base}_tagAH.vcf.gz"
        File output_vcf_idx = "~{vcf_base}_tagAH.vcf.gz.tbi"
    }

    runtime {
        docker: docker_image
        cpu: threads
        memory: memSizeGB + " GB"
        disks: "local-disk " + diskSizeGB + " SSD"
    }
}

task tag_haplotype_with_quality_filtering {
    input {
        File germline_VCF_input
        File germline_IDX_input
        File somatic_VCF_input
        File somatic_IDX_input
        File BAM
        File? BAI

        Int agreeing_gv_threshold
        Int disagreeing_gv_threshold
        
        String sample
        String vcf_base = basename(somatic_VCF_input, ".vcf.gz")


        String docker_image = "jiminpark/deepsomatic_postprocess:v2"
        Int threads = 20
        Int memSizeGB = 4
        Int diskSizeGB = round(size(BAM, 'G')) * 2
        # may need more disk than before..
    }

    command <<<
        if [[ "~{BAI}" == "" ]]
        then
                samtools index -@ ~{threads} ~{BAM}
        fi

        # echo each line of the script to stdout so we can see what is happening
        # to turn off echo do 'set +o xtrace'
        set -o xtrace

        python3 /opt/scripts/tag_haplotype_exclusive_multiprocessing.py -bam ~{BAM} -v ~{somatic_VCF_input} -g ~{germline_VCF_input} -d "output_files" -o "output_files/~{vcf_base}_tagAH.vcf.gz" -k ~{agreeing_gv_threshold} -m ~{disagreeing_gv_threshold} --max_workers ~{threads}
    >>>

    output {
        File output_vcf = "output_files/~{vcf_base}_tagAH.vcf.gz"
        File output_vcf_idx = "output_files/~{vcf_base}_tagAH.vcf.gz.tbi"
    }

    runtime {
        docker: docker_image
        cpu: threads
        memory: memSizeGB + " GB"
        disks: "local-disk " + diskSizeGB + " SSD"
    }
}