version 1.0

# https://github.com/nanoporegenomics/napu_wf/blob/a11bd9249d5cabf0feca2aa473f28786f75adbe7/wdl/tasks/dv-margin.wdl#L135
task merge_vcfs {
        input {
                Array[File] vcf_files
                Int memSizeGb = 64   # this can probably be reduced again..
                Int diskSizeGb = 5 * round(size(vcfFiles, 'G')) + 500
                String output_prefix
                String docker_image = "jiminpark/deepsomatic_postprocess"
        }

        command <<<
                set -o pipefail
                set -e
                set -u
                set -o xtrace

                mkdir bcftools.tmp

                # vcf merging
                bcftools concat -n ~{sep=" " vcf_files} | bcftools sort -T bcftools.tmp -O z -o ~{output_prefix}_merged.vcf.gz -
                bcftools index -t ~{output_prefix}_merged.vcf.gz.vcf.gz
        >>>

        output {
                File merged_vcf = "~{output_prefix}_merged.vcf.gz.vcf.gz"
                File merged_vcf_idx = "~{output_prefix}_merged.vcf.gz.vcf.gz.tbi"
        }

        runtime {
                docker: docker_image
                cpu: threads
                memory: memSizeGB + " GB"
                disks: "local-disk " + diskSizeGB + " SSD"
        }
}