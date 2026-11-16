version 1.0

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

                # restricts make_examples to this region (e.g. "chr20"); set it when tumor_bam is a per-region split
                String? region

                # extra run_deepsomatic flags, e.g. --use_default_pon_filtering=true
                String? additional_args

                String docker_image
                # 32 threads uses cores more fully than 64 (make_examples ran ~50-70% busy at 64)
                Int threads = 32
                # ~2 GB per thread, the same ratio as the old 64 threads / 128 GB
                Int memSizeGB = threads * 2
                Int diskSizeGB = 5 * round(size(tumor_bam, 'G')) + 50
                Int preemptible = 3
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

                run_deepsomatic \
                --model_type="~{model_type}" \
                --ref="~{reference}" \
                --reads_tumor=tumor.bam \
                --output_vcf="~{output_prefix}_tumor-only.vcf.gz" \
                --num_shards="~{threads}" \
                --logging_dir="log_outputs" \
                --sample_name_tumor="~{sample_name_tumor}" \
                ~{"--regions=" + region} \
                ${ADDITIONAL_ARGS}
        >>>

        output {
                File output_vcf = "~{output_prefix}_tumor-only.vcf.gz"
                File output_vcf_idx = "~{output_prefix}_tumor-only.vcf.gz.tbi"
        }

        runtime {
                preemptible: preemptible
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
                # region extraction is I/O-bound; more threads mostly sit idle during localization
                Int threads = 4
                Int memSizeGB = 8
                # full BAM + one region's slice (chr1 is ~8% of the genome)
                Int diskSizeGB = ceil(1.2 * size(bam, 'GB')) + 20
                Int preemptible = 3
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
                preemptible: preemptible
                docker: docker_image
                cpu: threads
                memory: memSizeGB + " GB"
                disks: "local-disk " + diskSizeGB + " SSD"
        }
}

# DRAFT: like split_bams, but reads only the region from gs:// instead of localizing the whole BAM.
# Requires a samtools built with GCS support; google/deepsomatic:1.9.0 has one (samtools 1.15.1, GCS=yes),
# jiminpark/deepsomatic_postprocess does not (GCS=no). Local (non-gs://) paths also work.
task split_bams_stream {
        input {
                File bam
                File bai
                String sample

                String region

                String docker_image = "google/deepsomatic:1.9.0"
                Int threads = 4
                Int memSizeGB = 8
                # only the region's slice is written locally (chr1 is ~8% of the genome)
                Int diskSizeGB = ceil(0.15 * size(bam, 'GB')) + 20
                Int preemptible = 3
        }

        parameter_meta {
                bam: { localization_optional: true }
                bai: { localization_optional: true }
        }

        command <<<
                set -euo pipefail
                set -o xtrace

                if [[ "~{bam}" == gs://* ]]
                then
                        # token from the VM's service account; keep it out of the xtrace log
                        set +o xtrace
                        GCS_OAUTH_TOKEN=$(curl -sf -H "Metadata-Flavor: Google" \
                                "http://metadata.google.internal/computeMetadata/v1/instance/service-accounts/default/token" \
                                | python3 -c 'import sys, json; print(json.load(sys.stdin)["access_token"])')
                        export GCS_OAUTH_TOKEN
                        set -o xtrace
                        export CURL_CA_BUNDLE=/etc/ssl/certs/ca-certificates.crt
                fi

                samtools view -bh -@ ~{threads} -X -o ~{sample}_~{region}.bam "~{bam}" "~{bai}" ~{region}
                samtools index -@ ~{threads} ~{sample}_~{region}.bam
        >>>

        output {
                File region_bam = "~{sample}_~{region}.bam"
                File region_bai = "~{sample}_~{region}.bam.bai"
        }

        runtime {
                preemptible: preemptible
                docker: docker_image
                cpu: threads
                memory: memSizeGB + " GB"
                disks: "local-disk " + diskSizeGB + " SSD"
        }
}