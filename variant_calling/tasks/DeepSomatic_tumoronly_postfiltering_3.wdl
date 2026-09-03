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

task tag_gnomad_filter_all {
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

        # keep all variant entries w/o gnomad tag
        bcftools view -i 'INFO/AF="."' ~{vcf_base}_gnomad.vcf.gz | bgzip  > ~{vcf_base}_gnomadfilterall.vcf.gz

        bcftools index -t  ~{vcf_base}_gnomadfilterall.vcf.gz

    >>>

    output {
        File gnomad_VCF = "~{vcf_base}_gnomad.vcf.gz"
        File gnomad_IDX = "~{vcf_base}_gnomad.vcf.gz.tbi"
        File gnomad_filter_VCF = "~{vcf_base}_gnomadfilterall.vcf.gz"
        File gnomad_filter_IDX = "~{vcf_base}_gnomadfilterall.vcf.gz.tbi"
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



task intersect_dv_and_gnomad {
    input {
        File germline_VCF
        File germline_IDX

        File gnomad_VCF
        File gnomad_IDX

        String sample
        String docker_image = "jiminpark/deepsomatic_postprocess:v3"

        Int threads = 1
        Int memSizeGB = 4
        Int diskSizeGB = 64
    }

    command <<<
        # echo each line of the script to stdout so we can see what is happening
        # to turn off echo do 'set +o xtrace'
        set -o xtrace

        mkdir bcftools_isec_temp
        bcftools isec -p bcftools_isec_temp ~{germline_VCF} ~{gnomad_VCF}
        cd bcftools_isec_temp

        mv 0002.vcf ~{sample}_dv_intersect_gnomad_0_01.vcf
        bgzip ~{sample}_dv_intersect_gnomad_0_01.vcf
        bcftools index -t ~{sample}_dv_intersect_gnomad_0_01.vcf.gz

        mv ~{sample}_dv_intersect_gnomad_0_01.vcf.gz ../
        mv ~{sample}_dv_intersect_gnomad_0_01.vcf.gz.tbi ../

    >>>

    output {
        File output_vcf = "~{sample}_dv_intersect_gnomad_0_01.vcf.gz"
        File output_vcf_idx = "~{sample}_dv_intersect_gnomad_0_01.vcf.gz.tbi"
    }

    runtime {
        docker: docker_image
        cpu: threads
        memory: memSizeGB + " GB"
        disks: "local-disk " + diskSizeGB + " SSD"
    }
}


task extract_snps_from_harmphase {
    input {
        File harmphase_VCF
        File? harmphase_IDX

        String sample
        String vcf_base = basename(harmphase_VCF, ".vcf.gz")

        String docker_image = "jiminpark/deepsomatic_postprocess"
        Int threads = 1
        Int memSizeGB = 4
        Int diskSizeGB = 64
    }

    command <<<
        set -o xtrace

        ~{if !defined(harmphase_IDX) then "bcftools index -t " + harmphase_VCF else ""}

        bcftools norm -m -any ~{harmphase_VCF} | bcftools view -v snps | bgzip > ~{vcf_base}_snps.vcf.gz
        bcftools index -t ~{vcf_base}_snps.vcf.gz
    >>>

    output {
        File output_vcf = "~{vcf_base}_snps.vcf.gz"
        File output_vcf_idx = "~{vcf_base}_snps.vcf.gz.tbi"
    }

    runtime {
        docker: docker_image
        cpu: threads
        memory: memSizeGB + " GB"
        disks: "local-disk " + diskSizeGB + " SSD"
    }
}

task tag_HQ {
    input {
        String script
        File bam
        File bai

        File germline_dv_intersect_gnomad
        File germline_dv_intersect_gnomad_idx
        File somatic_VCF_input
        File somatic_IDX_input

        Int agreeing_gv_threshold = 5
        Int disagreeing_gv_threshold = 0
        Int unphased_threshold = 10
        String? chromosome

        String vcf_base = basename(somatic_VCF_input, ".vcf.gz")
        
        String docker_image = "jiminpark/deepsomatic_postprocess:v3"
        Int threads = 30
        Int memSizeGB = 32
        Int diskSizeGB = round(size(bam, 'G')) * 4
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

        python3 /opt/scripts/~{script} \
        -bam ~{bam} \
        -v ~{somatic_VCF_input} \
        -g ~{germline_dv_intersect_gnomad} \
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

task filter_hap_exclusive {
    input {
        File somatic_VCF_input
        File somatic_IDX_input

        String sample
        String vcf_base = basename(somatic_VCF_input, ".vcf.gz")

        String docker_image = "jiminpark/deepsomatic_postprocess:v3"
        Int threads = 1
        Int memSizeGB = 4
        Int diskSizeGB = 64
    }

    command <<<
        set -o xtrace
        bcftools view -i 'INFO/AH=1' ~{somatic_VCF_input} | bgzip > ~{vcf_base}_AH1_only.vcf.gz
        bcftools index -t ~{vcf_base}_AH1_only.vcf.gz
    >>>

    output {
        File output_vcf = "~{vcf_base}_AH1_only.vcf.gz"
        File output_vcf_idx = "~{vcf_base}_AH1_only.vcf.gz.tbi"
    }

    runtime {
        docker: docker_image
        cpu: threads
        memory: memSizeGB + " GB"
        disks: "local-disk " + diskSizeGB + " SSD"
    }
}

task filter_hap_vaf {
    input {
        File somatic_VCF_input
        File somatic_IDX_input

        String sample
        String vcf_base = basename(somatic_VCF_input, ".vcf.gz")

        String docker_image = "jiminpark/deepsomatic_postprocess:v3"
        Int threads = 1
        Int memSizeGB = 4
        Int diskSizeGB = 64
    }

    command <<<
        set -o xtrace
        bcftools view -e '(INFO/HAP1_VAF_HQ>=0.8 && INFO/HAP2_VAF_HQ=0) || (INFO/HAP1_VAF_HQ=0 && INFO/HAP2_VAF_HQ>=0.8)' ~{somatic_VCF_input} | bgzip > ~{vcf_base}_hapVAF_0.8.vcf.gz
        bcftools index -t ~{vcf_base}_hapVAF_0.8.vcf.gz
    >>>

    output {
        File output_vcf = "~{vcf_base}_hapVAF_0.8.vcf.gz"
        File output_vcf_idx = "~{vcf_base}_hapVAF_0.8.vcf.gz.tbi"
    }

    runtime {
        docker: docker_image
        cpu: threads
        memory: memSizeGB + " GB"
        disks: "local-disk " + diskSizeGB + " SSD"
    }
}

# subtract tandem repeats and homopolymer regions
task subtract_trhp_regions {
    input {
        File somatic_VCF_input
        File somatic_IDX_input

        String sample
        String vcf_base = basename(somatic_VCF_input, ".vcf.gz")

        String docker_image = "jiminpark/deepsomatic_postprocess:v3"
        Int threads = 1
        Int memSizeGB = 4
        Int diskSizeGB = 64
    }

    command <<<
        # echo each line of the script to stdout so we can see what is happening
        # to turn off echo do 'set +o xtrace'
        set -o xtrace

        bedtools subtract -header -a ~{somatic_VCF_input} -b /opt/scripts/GRCh38_AllTandemRepeatsandHomopolymers_slop5.bed | bgzip > ~{vcf_base}_AllTandemRepeatsandHomopolymers.vcf.gz
        bcftools index -t ~{vcf_base}_AllTandemRepeatsandHomopolymers.vcf.gz
    >>>

    output {
        File output_vcf = "~{vcf_base}_AllTandemRepeatsandHomopolymers.vcf.gz"
        File output_vcf_idx = "~{vcf_base}_AllTandemRepeatsandHomopolymers.vcf.gz.tbi"
    }

    runtime {
        docker: docker_image
        cpu: threads
        memory: memSizeGB + " GB"
        disks: "local-disk " + diskSizeGB + " SSD"
    }
}