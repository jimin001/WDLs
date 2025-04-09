version 1.0

workflow DeepSomatic {
        input{
                String output_prefix
        }

        call deepSomatic {
                input:
                        output_prefix = output_prefix
        }

        call postProcess {
                input:
                        vcf = deepSomatic.output_vcf,
                        output_prefix = output_prefix
        }
}

################### tasks #####################

task deepSomatic {
        input {

                File reference
                File reference_idx
                File normal_bam
                File normal_bam_idx
                File tumor_bam
                File tumor_bam_idx

                String output_prefix
                String log_dir_path
                String sample_name_tumor
                String sample_name_normal
                String docker_image
                String model_type

                # --regions=chr20
                String? additional_args

                # must be a tar.gz of FILES, not directory
                # File? model_file_tar

                File? model_file
                File? model_file_idx
                File? model_file_example

                # example: "weights-422-0.976350.ckpt"
                String? custom_model = sub(split(model_file_idx, "\\.")[0], "/$", "")


                Int memSizeGB = 128
                Int threadCount = 64
                Int diskSizeGB = 128

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


                if [[ "~{model_file}" == "" ]]
                then
                        run_deepsomatic \
                        --model_type="~{model_type}" \
                        --ref="~{reference}" \
                        --reads_normal="~{normal_bam}" \
                        --reads_tumor="~{tumor_bam}" \
                        --output_vcf="~{output_prefix}.vcf.gz" \
                        --num_shards=~{threadCount} \
                        --logging_dir="~{log_dir_path}" \
                        --sample_name_tumor="~{sample_name_tumor}" \
                        --sample_name_normal="~{sample_name_normal}" \
                        ${ADDITIONAL_ARGS}
                else
                        run_deepsomatic \
                        --model_type="~{model_type}" \
                        --ref="~{reference}" \
                        --reads_normal="~{normal_bam}" \
                        --reads_tumor="~{tumor_bam}" \
                        --output_vcf="~{output_prefix}.vcf.gz" \
                        --num_shards=~{threadCount} \
                        --logging_dir="~{log_dir_path}" \
                        --sample_name_tumor="~{sample_name_tumor}" \
                        --sample_name_normal="~{sample_name_normal}" \
                        --customized_model=~{custom_model} \
                        ${ADDITIONAL_ARGS}
                fi
        >>>

        output {
                File output_vcf = "~{output_prefix}.vcf.gz"
        }

        runtime {
                memory: memSizeGB + " GB"
                cpu: threadCount
                disks: "local-disk " + diskSizeGB + " SSD"
                docker: docker_image
                preemptible: 1
        }
}


task postProcess {
        input {
                File vcf
                String output_prefix

                Int memSizeGB = 4
                Int diskSizeGB = 128
                Int threadCount = 2

                String docker_image = "jiminpark/tools:latest"
        }

        command <<<
                bcftools filter -i 'FILTER="PASS"' "~{vcf}" | bgzip > "~{output_prefix}.somatic_only.vcf.gz"

                tabix -p vcf "~{output_prefix}.somatic_only.vcf.gz"
        >>>

        output{
                File somatic_only_vcf = "~{output_prefix}.somatic_only.vcf.gz"
                File somatic_only_vcf_tbi = "~{output_prefix}.somatic_only.vcf.gz.tbi"
        }

        runtime {
                memory: memSizeGB + " GB"
                cpu: threadCount
                disks: "local-disk " + diskSizeGB + " SSD"
                docker: docker_image
                preemptible: 1
        }

}

