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

task filter_out_germline_variants {
    input {
        File germline_VCF_pass_only
        File germline_IDX_pass_only
        File somatic_VCF_pass_only
        File somatic_IDX_pass_only

        String sample

        String docker_image = "jiminpark/deepsomatic_postprocess"
        Int threads = 1
        Int memSizeGB = 4
        Int diskSizeGB = 64
    }

    command <<<
        # echo each line of the script to stdout so we can see what is happening
        # to turn off echo do 'set +o xtrace'
        set -o xtrace

        SAMPLE="~{sample}"
        mkdir ${SAMPLE}
        bcftools isec -p ${SAMPLE} ~{somatic_VCF_pass_only} ~{germline_VCF_pass_only}

        cd ${SAMPLE}
        mv 0000.vcf ~{sample}_deepsomatic_only.vcf
        mv 0001.vcf ~{sample}_deepvariant_only.vcf
        mv 0002.vcf ~{sample}_intersection.vcf
        rm 0003.vcf
        rm *txt

        bgzip ~{sample}_deepsomatic_only.vcf
        bcftools index -t ~{sample}_deepsomatic_only.vcf.gz

        mv ~{sample}_deepsomatic_only.vcf.gz ../
        mv ~{sample}_deepsomatic_only.vcf.gz.tbi ../

    >>>

    output {
        File deepsomatic_only_vcf = "~{sample}_deepsomatic_only.vcf.gz"
        File deepsomatic_only_vcf_idx = "~{sample}_deepsomatic_only.vcf.gz.tbi"
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
        File deepsomatic_only_VCF
        File deepsomatic_only_IDX

        String sample

        String docker_image = "jiminpark/deepsomatic_postprocess"
        Int threads = 1
        Int memSizeGB = 4
        Int diskSizeGB = 64
    }

    command <<<
        # echo each line of the script to stdout so we can see what is happening
        # to turn off echo do 'set +o xtrace'
        set -o xtrace

        bcftools view -e 'FMT/GQ<20' ~{deepsomatic_only_VCF} | bcftools view -e 'FMT/DP<10' | bgzip > ~{sample}_deepsomatic_only_GQ20_DP10.vcf.gz
        bcftools index -t ~{sample}_deepsomatic_only_GQ20_DP10.vcf.gz
    >>>

    output {
        File deepsomatic_only_GQ20_DP10_vcf = "~{sample}_deepsomatic_only_GQ20_DP10.vcf.gz"
        File deepsomatic_only_GQ20_DP10_vcf_idx = "~{sample}_deepsomatic_only_GQ20_DP10.vcf.gz.tbi"
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
        File deepsomatic_only_GQ20_DP10_VCF
        File deepsomatic_only_GQ20_DP10_IDX

        String sample

        String docker_image = "jiminpark/deepsomatic_postprocess"
        Int threads = 1
        Int memSizeGB = 4
        Int diskSizeGB = 64
    }

    command <<<
        # echo each line of the script to stdout so we can see what is happening
        # to turn off echo do 'set +o xtrace'
        set -o xtrace

        # bcftools view -h ~{deepsomatic_only_GQ20_DP10_VCF} > ~{sample}_header_temp.vcf
        bedtools subtract -header -a ~{deepsomatic_only_GQ20_DP10_VCF} -b /opt/scripts/GRCh38_segdups.bed | bgzip > ~{sample}_deepsomatic_only_GQ20_DP10_segdup.vcf.gz
        bcftools index -t ~{sample}_deepsomatic_only_GQ20_DP10_segdup.vcf.gz
    >>>

    output {
        File deepsomatic_only_GQ20_DP10_segdup_vcf = "~{sample}_deepsomatic_only_GQ20_DP10_segdup.vcf.gz"
        File deepsomatic_only_GQ20_DP10_segdup_idx = "~{sample}_deepsomatic_only_GQ20_DP10_segdup.vcf.gz.tbi"
    }

    runtime {
        docker: docker_image
        cpu: threads
        memory: memSizeGB + " GB"
        disks: "local-disk " + diskSizeGB + " SSD"
    }
}

task tag_haplotype {
    input {
        File deepsomatic_only_GQ20_DP10_segdup_VCF
        File deepsomatic_only_GQ20_DP10_segdup_IDX
        File BAM
        File BAI
        
        String sample

        String docker_image = "jiminpark/deepsomatic_postprocess"
        Int threads = 1
        Int memSizeGB = 4
        Int diskSizeGB = 128
    }

    command <<<
        # echo each line of the script to stdout so we can see what is happening
        # to turn off echo do 'set +o xtrace'
        set -o xtrace

        python3 /opt/scripts/haplotypeonly_somatic_vcf.py -v ~{deepsomatic_only_GQ20_DP10_segdup_VCF} -bam ~{BAM} -o ~{sample}_deepsomatic_only_GQ20_DP10_segdup_tagAH.vcf -u 10
        bgzip ~{sample}_deepsomatic_only_GQ20_DP10_segdup_tagAH.vcf
        bcftools index -t ~{sample}_deepsomatic_only_GQ20_DP10_segdup_tagAH.vcf.gz
    >>>

    output {
        File deepsomatic_only_GQ20_DP10_segdup_tagAH_vcf = "~{sample}_deepsomatic_only_GQ20_DP10_segdup_tagAH.vcf.gz"
        File deepsomatic_only_GQ20_DP10_segdup_tagAH_idx = "~{sample}_deepsomatic_only_GQ20_DP10_segdup_tagAH.vcf.gz.tbi"
    }

    runtime {
        docker: docker_image
        cpu: threads
        memory: memSizeGB + " GB"
        disks: "local-disk " + diskSizeGB + " SSD"
    }
}

task tag_gnomad {
    input {
        File deepsomatic_only_GQ20_DP10_segdup_tagAH_VCF
        File deepsomatic_only_GQ20_DP10_segdup_tagAH_IDX

        String sample

        String docker_image = "jiminpark/snpeff"
        Int threads = 1
        Int memSizeGB = 4
        Int diskSizeGB = 128
    }

    command <<<
        # echo each line of the script to stdout so we can see what is happening
        # to turn off echo do 'set +o xtrace'
        set -o xtrace

        java -jar /home/apps/snpEff/SnpSift.jar Annotate -noId /opt/scripts/gnomad.genomes.v4.1.sites.small.vcf.bgz \
        ~{deepsomatic_only_GQ20_DP10_segdup_tagAH_VCF} | bgzip > ~{sample}_deepsomatic_only_GQ20_DP10_segdup_tagAH_gnomad.vcf.gz

        bcftools index -t ~{sample}_deepsomatic_only_GQ20_DP10_segdup_tagAH_gnomad.vcf.gz

        bcftools view -i 'INFO/AF<0.001 | INFO/AF="."' ~{sample}_deepsomatic_only_GQ20_DP10_segdup_tagAH_gnomad.vcf.gz \
        | bcftools view -i 'INFO/AH=1' | bgzip  > ~{sample}_deepsomatic_only_GQ20_DP10_segdup_tagAH_gnomadrare_AH1.vcf.gz

        bcftools index -t ~{sample}_deepsomatic_only_GQ20_DP10_segdup_tagAH_gnomadrare_AH1.vcf.gz

    >>>

    output {
        File gnomad_VCF = "~{sample}_deepsomatic_only_GQ20_DP10_segdup_tagAH_gnomad.vcf.gz"
        File gnomad_IDX = "~{sample}_deepsomatic_only_GQ20_DP10_segdup_tagAH_gnomad.vcf.gz.tbi"
        File gnomad_filter_VCF = "~{sample}_deepsomatic_only_GQ20_DP10_segdup_tagAH_gnomadrare_AH1.vcf.gz"
        File gnomad_filter_IDX = "~{sample}_deepsomatic_only_GQ20_DP10_segdup_tagAH_gnomadrare_AH1.vcf.gz.tbi"
    }

    runtime {
        docker: docker_image
        cpu: threads
        memory: memSizeGB + " GB"
        disks: "local-disk " + diskSizeGB + " SSD"
    }
}