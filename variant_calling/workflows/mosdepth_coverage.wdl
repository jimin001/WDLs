version 1.0

# Binned read depth for one BAM with mosdepth, e.g. for coverage outlines under call strip plots.
# Run per sample (one Terra data-table row). --fast-mode skips per-base CIGAR/mate-overlap
# handling, which only matters at base resolution, and --no-per-base drops the per-base output,
# so the result is mean depth per bin_size window plus whole-genome summaries.

workflow mosdepthCoverage {
    input {
        String sample_name
        File bam
        File bam_index

        # 1 Mb matches the ONT perbase outline bins in the strip plot
        Int bin_size = 1000000
        # 20 matches check_pacbio_signal and DeepSomaticFilterNN on Illumina; drops multi-mapped
        # reads so low-mappability regions show as dips. 0 counts every mapped read
        Int min_mapq = 20
        # unmapped, secondary, QC-fail, duplicate and supplementary; mosdepth's default 1796 keeps
        # supplementary, while ONT perbase (-F 3848) drops it
        Int exclude_flags = 3844
    }

    call mosdepth {
        input:
            sample_name = sample_name,
            bam = bam,
            bam_index = bam_index,
            bin_size = bin_size,
            min_mapq = min_mapq,
            exclude_flags = exclude_flags
    }

    output {
        File regions_bed = mosdepth.regions_bed
        File regions_bed_index = mosdepth.regions_bed_index
        File summary = mosdepth.summary
        File global_dist = mosdepth.global_dist
    }
}

task mosdepth {
    input {
        String docker_image = "quay.io/biocontainers/mosdepth:0.3.12--h0ec343a_0"

        String sample_name
        File bam
        File bam_index

        Int bin_size = 1000000
        Int min_mapq = 20
        Int exclude_flags = 3844

        # mosdepth uses threads only for BAM decompression; gains flatten past ~4
        Int threads = 4
        Int memSizeGB = 8
        # only the BAM is localized and nothing large is written, so a restart loses little
        Int preemptible = 2
    }
    Int diskSizeGB = ceil(size(bam, "GB")) + 20

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

        # index must sit next to the BAM; Cromwell may localize them to different dirs
        ln -s ~{bam} input.bam
        ln -s ~{bam_index} input.bam.bai

        mosdepth \
            --threads ~{threads} \
            --fast-mode \
            --no-per-base \
            --by ~{bin_size} \
            --mapq ~{min_mapq} \
            --flag ~{exclude_flags} \
            ~{sample_name} \
            input.bam

        cat ~{sample_name}.mosdepth.summary.txt
    >>>

    output {
        File regions_bed = "working/~{sample_name}.regions.bed.gz"
        File regions_bed_index = "working/~{sample_name}.regions.bed.gz.csi"
        File summary = "working/~{sample_name}.mosdepth.summary.txt"
        File global_dist = "working/~{sample_name}.mosdepth.global.dist.txt"
    }

    runtime {
        memory: memSizeGB + " GB"
        cpu: threads
        # HDD: the BAM is read once sequentially, so SSD adds cost without much speedup
        disks: "local-disk " + diskSizeGB + " HDD"
        docker: docker_image
        preemptible: preemptible
        # retry transient failures such as registry pull errors
        maxRetries: 2
    }
}
