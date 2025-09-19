version 1.0
import "merge_vcfs.wdl" as merge_vcfs
import "DeepSomatic_tumoronly_postfiltering.wdl" as filter


workflow DeepSomatic {
        input{
            Array[File] bam_files
            Array[File] bam_files_idx
            String output_prefix
        }

        scatter (bam_file in zip(bam_files, bam_files_idx)) {
            call deepSomatic {
                input:
                        tumor_bam = bam_file.left,
                        tumor_bam_idx = bam_file.right,
                        output_prefix = output_prefix
            }
        }

        call merge_vcfs.merge_vcfs as merge_vcfs {
                input:
                        vcf_files = deepSomatic.output_vcf,
                        output_prefix = output_prefix
        }        

        output {
            File DeepSomaticVCF = merge_vcfs.merged_vcf
            File DeepSomaticVCFIndex = merge_vcfs.merged_vcf_idx
        }
}

################### tasks #####################

task deepSomatic {
        input {

                File reference
                File reference_idx
                File tumor_bam
                File tumor_bam_idx

                String output_prefix
                String sample_name_tumor

                # "WGS" or "ONT_R104"
                String model_type

                # --regions=chr20
                String? additional_args

                String docker_image
                Int threads
                Int memSizeGB = 128
                Int diskSizeGB = 256
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

                mkdir log_outputs

                run_deepsomatic \
                --model_type="~{model_type}" \
                --ref="~{reference}" \
                --reads_tumor="~{tumor_bam}" \
                --output_vcf="~{output_prefix}_tumor-only.vcf.gz" \
                --num_shards="~{threads}" \
                --logging_dir="log_outputs" \
                --sample_name_tumor="~{sample_name_tumor}" \
                ${ADDITIONAL_ARGS}
        >>>

        output {
                File output_vcf = "~{output_prefix}_tumor-only.vcf.gz"
                File output_vcf_idx = "~{output_prefix}_tumor-only.vcf.gz.tbi"
        }

        runtime {
                preemtiple: 2
                docker: docker_image
                cpu: threads
                memory: memSizeGB + " GB"
                disks: "local-disk " + diskSizeGB + " SSD"
        }
}

