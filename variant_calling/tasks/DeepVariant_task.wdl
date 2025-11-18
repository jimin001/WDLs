version 1.0

task deepVariant {
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
                Int diskSizeGB = 5 * round(size(tumor_bam, 'G')) + 50
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

                ln -s ~{tumor_bam} tumor.bam
                ln -s ~{tumor_bam_idx} tumor.bam.bai

                if [[ "~{additional_args}" == "" ]]
                then
                        ADDITIONAL_ARGS=""
                else
                        ADDITIONAL_ARGS="~{additional_args}"
                fi

                mkdir log_outputs

                run_deepvariant \
                --model_type="~{model_type}" \
                --ref="~{reference}" \
                --reads=tumor.bam \
                --output_vcf="~{output_prefix}_deepvariant.vcf.gz" \
                --num_shards="~{threads}" \
                --logging_dir="log_outputs" \
                ${ADDITIONAL_ARGS}
        >>>

        output {
                File output_vcf = "~{output_prefix}_deepvariant.vcf.gz"
                File output_vcf_idx = "~{output_prefix}_deepvariant.vcf.gz.tbi"
        }

        runtime {
                preemptible: 2
                docker: docker_image
                cpu: threads
                memory: memSizeGB + " GB"
                disks: "local-disk " + diskSizeGB + " SSD"
        }
}

task split_bams {
        input {
                File bam
                File? bai
                String sample

                String region

                String docker_image = "jiminpark/deepsomatic_postprocess"
                Int threads
                Int memSizeGB = 128
                Int diskSizeGB = 2 * round(size(bam, 'G')) + 50
        }

        command <<<
                if [[ "~{bai}" == "" ]]
                then
                        samtools index -@ ~{threads} ~{bam}
                else
                        continue
                fi

                samtools view -bh -@ ~{threads} ~{bam} ~{region} > ~{sample}_~{region}.bam
                samtools index -@ ~{threads} ~{sample}_~{region}.bam
        >>>

        output {
                File region_bam = "~{sample}_~{region}.bam"
                File region_bai = "~{sample}_~{region}.bam.bai"
        }

        runtime {
                preemptible: 2
                docker: docker_image
                cpu: threads
                memory: memSizeGB + " GB"
                disks: "local-disk " + diskSizeGB + " SSD"
        }
}