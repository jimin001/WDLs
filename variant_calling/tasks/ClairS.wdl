version 1.0

workflow ClairS {
        input{
                File reference
                File reference_idx
                String output_prefix
        }

        call ClairS {
                input:
                        reference = reference,
                        reference_idx = reference_idx
        }

        call postProcessClairS {
                input:
                        output_prefix = output_prefix,
                        output_vcf = ClairS.output_vcf,
                        output_vcf_idx = ClairS.output_vcf_idx,
                        indel_vcf = ClairS.indel_vcf,
                        indel_vcf_idx = ClairS.indel_vcf_idx
        }

        #call sompy {
        #        input:
        #                reference = reference,
        #                reference_idx = reference_idx,
        #                snv_indel_vcf = postProcessClairS.somatic_only_vcf,
        #                snv_indel_vcf_idx = postProcessClairS.somatic_only_vcf_tbi
        #}
}

################### tasks #####################

task ClairS {
        input {

                File reference
                File reference_idx
                File normal_bam
                File normal_bam_idx
                File tumor_bam
                File tumor_bam_idx

                String sample_name_tumor
                String sample_name_normal

                Int threads = 64

                # "WGS" or "ONT_R104"
                String platform

                # for HiFi and ONT: --enable_indel_calling
                # --region chr1
                String? additional_args

                String docker_image = "hkubal/clairs:latest"
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

                if [[ "~{additional_args}" == "" ]]
                then
                        ADDITIONAL_ARGS=""
                else
                        ADDITIONAL_ARGS="~{additional_args}"
                fi

                # mkdir clairs_outputs
                # cd clairs_outputs

                /opt/bin/run_clairs \
                --tumor_bam_fn="~{tumor_bam}" \
                --normal_bam_fn="~{normal_bam}" \
                --ref_fn="~{reference}" \
                --threads="~{threads}" \
                --platform="~{platform}" \
                --output_dir="clairs_outputs" \
                --sample_name "ClairS" \
                ${ADDITIONAL_ARGS}
        >>>

        output {
                File output_vcf = "clairs_outputs/output.vcf.gz"
                File output_vcf_idx = "clairs_outputs/output.vcf.gz.tbi"
                File indel_vcf = "clairs_outputs/indel.vcf.gz"
                File indel_vcf_idx = "clairs_outputs/indel.vcf.gz.tbi"
        }

        runtime {
                docker: docker_image
        }
}

task postProcessClairS {
        input {
                String output_prefix

                # files output from task "ClairS"
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

task sompy {
        input {
                String docker_image = "pkrusche/hap.py:latest"

                File truth_vcf
                File truth_vcf_idx

                File snv_indel_vcf
                File snv_indel_vcf_idx

                File reference
                File reference_idx
                File bed_file

                String output_prefix

                # -l chr20
                String? additional_args

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

                if [[ "~{additional_args}" == "" ]]
                then
                        ADDITIONAL_ARGS=""
                else
                        ADDITIONAL_ARGS="~{additional_args}"
                fi

                # run sompy on snvs
                /opt/hap.py/bin/som.py -N \
                ~{truth_vcf} \
                ~{snv_indel_vcf} \
                --restrict-regions ~{bed_file} \
                -r ~{reference} \
                -o "~{output_prefix}_snv.indel.sompy" \
                --feature-table generic \
                ${ADDITIONAL_ARGS}


        >>>

        output {
                File features_snvs = "~{output_prefix}_snv.indel.sompy.features.csv"
                File metrics_snvs = "~{output_prefix}_snv.indel.sompy.metrics.json"
                File stats_snvs = "~{output_prefix}_snv.indel.sompy.stats.csv"
        }

        runtime {
                docker: docker_image
        }

}
