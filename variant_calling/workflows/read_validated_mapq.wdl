version 1.0

# Per-site mapping-quality summary, in matched Illumina and ONT BAMs, for a list of SNV sites
# (e.g. illumina_read_validated_ont_PoN_SNVs_illumina_call_status.tsv from
# scripts/plots/rush_illumina_read_validated_not_called.py).
#
# mapqAtSites runs once per BAM: GATK PrintReads -L streams only the reads overlapping the sites from GCS
# (the full BAM is never localized) into a subset BAM on the VM's local disk, then
# scripts/plots/rush_read_validated_mapq.py summarizes it. Only the per-BAM TSV is a task output, so the
# subset BAM is never copied to the workspace bucket. concatTsvs merges the per-BAM TSVs.
#
# GATK is used for the streaming because pysam in jiminpark/deepsomatic_postprocess:v3 cannot do gs:// over
# HTTPS (its bundled NSS libcurl fails with CA error 77) and that image's samtools 1.17 has no gs:// support.
# GATK reads gs:// with the VM's own credentials, so no token handling is needed. The GATK image's
# python 3.6 / pysam 0.15.3 runs the summary script.

workflow readValidatedMapq {
    input {
        # matched arrays, one entry per sample (e.g. this.samples.<column> on a Terra sample_set)
        Array[String] sample_names
        Array[File] illumina_bams
        Array[File] illumina_bam_indexes
        Array[File] ont_bams
        Array[File] ont_bam_indexes

        # TSV with header; needs columns sample, chrom, pos, ref, alt, illumina_status
        # local copy: /private/groups/patenlab/jimin/data/VCF/CARD/RUSH_Illumina_analysis/illumina_read_validated_ont_PoN_SNVs_illumina_call_status.tsv
        File sites_tsv

        # /private/groups/patenlab/jimin/scripts/plots/rush_read_validated_mapq.py
        File script

        # matches check_pacbio_signal; there is deliberately no MAPQ minimum
        Int min_baseq = 10
    }

    scatter (i in range(length(sample_names))) {
        call mapqAtSites as mapqIllumina {
            input:
                sample_name = sample_names[i],
                platform = "illumina",
                bam = illumina_bams[i],
                bam_index = illumina_bam_indexes[i],
                sites_tsv = sites_tsv,
                script = script,
                min_baseq = min_baseq
        }

        call mapqAtSites as mapqONT {
            input:
                sample_name = sample_names[i],
                platform = "ont",
                bam = ont_bams[i],
                bam_index = ont_bam_indexes[i],
                sites_tsv = sites_tsv,
                script = script,
                min_baseq = min_baseq
        }
    }

    call concatTsvs {
        input:
            tsvs = flatten([mapqIllumina.mapq_tsv, mapqONT.mapq_tsv])
    }

    output {
        File mapq_tsv = concatTsvs.merged_tsv
    }
}

task mapqAtSites {
    input {
        String docker_image = "us.gcr.io/broad-gatk/gatk:4.5.0.0"

        String sample_name
        # illumina or ont
        String platform

        File bam
        File bam_index
        File sites_tsv
        File script

        Int min_baseq = 10

        Int memSizeGB = 4
        Int diskSizeGB = 20
        Int preemptible = 2
    }

    parameter_meta {
        # passed through as gs:// URLs and streamed by GATK instead of being copied to the VM
        bam: { localization_optional: true }
        bam_index: { localization_optional: true }
    }

    String prefix = "~{sample_name}_~{platform}"

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

        mkdir working
        cd working

        # one chrom:pos-pos interval per site of this sample, located by header name
        awk -F'\t' -v s=~{sample_name} '
            NR == 1 { for (i = 1; i <= NF; i++) col[$i] = i; next }
            $col["sample"] == s { print $col["chrom"] ":" $col["pos"] "-" $col["pos"] }
        ' ~{sites_tsv} | sort -u > sites.intervals

        if [ ! -s sites.intervals ]
        then
            echo "ERROR: no sites for sample ~{sample_name} in ~{sites_tsv}" >&2
            exit 1
        fi
        wc -l sites.intervals

        # no read filters at all, so MAPQ-0, duplicate and secondary reads all reach the pileup
        # (which applies the same filters as check_pacbio_signal, minus the MAPQ minimum)
        gatk --java-options "-Xmx$(( ~{memSizeGB} * 1024 - 768 ))m" PrintReads \
            -I ~{bam} \
            --read-index ~{bam_index} \
            -L sites.intervals \
            --disable-tool-default-read-filters true \
            --create-output-bam-index true \
            -O subset.bam

        # bai column left empty: htslib finds subset.bai itself, and pysam 0.15.3 (GATK image) raises
        # "TypeError: expected bytes, str found" on fetch when index_filename is passed explicitly
        printf "sample\tplatform\tbam\tbai\n~{sample_name}\t~{platform}\tsubset.bam\t\n" > manifest.tsv

        python3 ~{script} \
            --sites ~{sites_tsv} \
            --manifest manifest.tsv \
            --out ~{prefix}.mapq.tsv \
            --min-baseq ~{min_baseq}
    >>>

    output {
        # subset.bam is deliberately not an output, so it never reaches the workspace bucket
        File mapq_tsv = "working/~{prefix}.mapq.tsv"
    }

    runtime {
        memory: memSizeGB + " GB"
        cpu: 1
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: docker_image
        preemptible: preemptible
        maxRetries: 2
    }
}

task concatTsvs {
    input {
        String docker_image = "ubuntu:22.04"

        Array[File] tsvs

        Int preemptible = 2
    }

    command <<<
        set -o pipefail
        set -e
        set -u
        set -o xtrace

        # header from the first file, body from all of them
        awk 'FNR == 1 && NR != 1 { next } { print }' ~{sep=" " tsvs} > read_validated_mapq.tsv
        wc -l read_validated_mapq.tsv
    >>>

    output {
        File merged_tsv = "read_validated_mapq.tsv"
    }

    runtime {
        memory: "1 GB"
        cpu: 1
        disks: "local-disk 10 SSD"
        docker: docker_image
        preemptible: preemptible
        # retry transient failures such as Docker Hub pull rate limits
        maxRetries: 2
    }
}
