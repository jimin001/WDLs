version 1.0

workflow Strelka2 {
        input{
                File reference
                File reference_idx
                String output_prefix
        }

        call Strelka2 {
                input:
                        reference = reference,
                        reference_idx = reference_idx
        }

        call postProcess {
                input:
                        output_prefix = output_prefix,
                        output_vcf = Strelka2.output_vcf,
                        output_vcf_idx = Strelka2.output_vcf_idx,
                        indel_vcf = Strelka2.indel_vcf,
                        indel_vcf_idx = Strelka2.indel_vcf_idx
        }
}

################### tasks #####################

task Strelka2 {
        input {

                File reference
                File reference_idx
                File normal_bam
                File normal_bam_idx
                File tumor_bam
                File tumor_bam_idx

                Int threads = 64
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

                STRELKA_INSTALL_PATH=/private/home/jpark621/software/strelka-2.9.10.centos6_x86_64/bin
                BED=/private/groups/patenlab/jimin/data/BED/strelka2_whole_genome_regions.bed.gz

                ${STRELKA_INSTALL_PATH}/configureStrelkaSomaticWorkflow.py \
                --normalBam "~{normal_bam}" \
                --tumorBam "~{tumor_bam}" \
                --referenceFasta "~{reference}" \
                --runDir "strelka_analysis_directory" \
                --callRegions ${BED}

                # configure step outputs "runWorkflow.py"
                ${STRELKA_ANALYSIS_PATH}/runWorkflow.py -m local -j 64
        >>>

        output {
                File output_vcf = "strelka_analysis_directory/results/variants/somatic.snvs.vcf.gz"
                File output_vcf_idx = "sstrelka_analysis_directory/results/variants/somatic.snvs.vcf.gz.tbi"
                File indel_vcf = "strelka_analysis_directory/results/variants/somatic.indels.vcf.gz"
                File indel_vcf_idx = "strelka_analysis_directory/results/variants/somatic.indels.vcf.gz.tbi"
        }
}

task postProcess {
        input {
                String output_prefix

                # files output from task "Strelka2"
                File output_vcf
                File output_vcf_idx
                File indel_vcf
                File indel_vcf_idx

                Int memSizeGB = 4
                Int diskSizeGB = 128
                Int threadCount = 2

                String docker_image = "jiminpark/tools:latest"
        }

        command <<<
                bcftools merge --force-samples "~{output_vcf}" "~{indel_vcf}" | bgzip > "~{output_prefix}.snv.indel.vcf.gz"
                bcftools index -t "~{output_prefix}.snv.indel.vcf.gz"

                bcftools filter -i 'FILTER="PASS"' "~{output_prefix}.snv.indel.vcf.gz" | bgzip > "~{output_prefix}.snv.indel.somatic_only.vcf.gz"

                bcftools index -t "~{output_prefix}.snv.indel.somatic_only.vcf.gz"
        >>>

        output{
                File somatic_only_vcf = "~{output_prefix}.snv.indel.somatic_only.vcf.gz"
                File somatic_only_vcf_tbi = "~{output_prefix}.snv.indel.somatic_only.vcf.gz.tbi"
        }

        runtime {
                memory: memSizeGB + " GB"
                cpu: threadCount
                disks: "local-disk " + diskSizeGB + " SSD"
                docker: docker_image
                preemptible: 1
        }

}
