version 1.0

# Sample swap check between matched ONT and Illumina data with somalier.
# extract runs on each ONT and Illumina BAM, relate runs on all digests together,
# then every expected ONT/Illumina pair is checked for a same-individual match.
#
# Use BAMs for both platforms, not variants-only VCFs: a VCF with no 0/0 calls
# leaves hom-ref sites unknown and inflates relatedness between unrelated samples (~0.67 seen on
# MOSR_006 vs MOSR_041 ONT VCFs).

workflow somalierSampleSwap {
    input {
        # matched arrays, one entry per sample (e.g. this.samples.<column> on a Terra sample_set)
        Array[String] sample_names
        Array[File] ont_bams
        Array[File] ont_bam_indexes
        Array[File] illumina_bams
        Array[File] illumina_bam_indexes

        # somalier sites with chr-prefixed contigs; the nochr version finds 0 sites on chr-named data
        # local copy: /private/groups/patenlab/jimin/data/VCF/CARD/resources/sites.hg38.chr.vcf.gz
        File sites_vcf

        # reference both platforms were aligned to (GRCh38 no_alt analysis set)
        File reference
        File reference_fai

        # same individual is ~1.0 with BAMs on both sides, unrelated ~0
        Float min_match_relatedness = 0.8
        Int min_shared_sites = 1000
    }

    scatter (i in range(length(sample_names))) {
        call somalierExtract as extractONT {
            input:
                sample_name = sample_names[i],
                platform = "ONT",
                bam = ont_bams[i],
                bam_index = ont_bam_indexes[i],
                reference = reference,
                reference_fai = reference_fai,
                sites_vcf = sites_vcf
        }

        call somalierExtract as extractIllumina {
            input:
                sample_name = sample_names[i],
                platform = "ILMN",
                bam = illumina_bams[i],
                bam_index = illumina_bam_indexes[i],
                reference = reference,
                reference_fai = reference_fai,
                sites_vcf = sites_vcf
        }
    }

    call somalierRelate {
        input:
            digests = flatten([extractONT.digest, extractIllumina.digest]),
            sample_names = sample_names,
            min_match_relatedness = min_match_relatedness,
            min_shared_sites = min_shared_sites
    }

    output {
        Array[File] ont_digests = extractONT.digest
        Array[File] illumina_digests = extractIllumina.digest
        File relate_pairs = somalierRelate.pairs
        File relate_samples = somalierRelate.samples
        File relate_html = somalierRelate.html
        File swap_report = somalierRelate.swap_report
        File unexpected_matches = somalierRelate.unexpected_matches
        String swap_check_status = somalierRelate.status
    }
}

task somalierExtract {
    input {
        String docker_image = "brentp/somalier:v0.3.5"

        String sample_name
        # ONT or ILMN
        String platform

        File bam
        File bam_index
        File reference
        File reference_fai
        File sites_vcf

        Int threads = 2
        Int memSizeGB = 4
    }
    Int file_size = ceil(size(bam, "GB") + size(reference, "GB"))
    Int diskSizeGB = file_size + 20

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

        # index and fai must sit next to their files; Cromwell may localize them to different dirs
        ln -s ~{bam} input.bam
        ln -s ~{bam_index} input.bam.bai
        ln -s ~{reference} ref.fa
        ln -s ~{reference_fai} ref.fa.fai

        # name the digest <platform>_<sample> from the inputs instead of the @RG SM tag:
        # some ONT BAMs have no read group, and SOMALIER_SAMPLE_NAME takes precedence over @RG
        export SOMALIER_SAMPLE_NAME=~{platform}_~{sample_name}
        somalier extract -s ~{sites_vcf} -f ref.fa -d . input.bam 2>&1 | tee extract.log

        if grep -q "found 0 sites" extract.log
        then
            echo "ERROR: 0 sites found; check that sites_vcf contigs match the reference (chr vs nochr)" >&2
            exit 1
        fi
    >>>

    output {
        File digest = "working/~{platform}_~{sample_name}.somalier"
    }

    runtime {
        memory: memSizeGB + " GB"
        cpu: threads
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: docker_image
        # retry transient failures such as Docker Hub pull rate limits
        maxRetries: 2
    }
}

task somalierRelate {
    input {
        String docker_image = "brentp/somalier:v0.3.5"

        Array[File] digests
        Array[String] sample_names

        Float min_match_relatedness = 0.8
        Int min_shared_sites = 1000

        Int threads = 2
        Int memSizeGB = 4
    }
    Int diskSizeGB = 20

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

        somalier relate -o somalier ~{sep=" " digests}

        # Status comes only from each sample's own ONT/Illumina pair. The best Illumina match and
        # unexpected_matches (different samples at >= threshold) are informational, to show where
        # a mismatched sample's data actually belongs. Digests are named <ONT|ILMN>_<sample>.
        # pairs.tsv columns: $1 sample_a, $2 sample_b, $3 relatedness, $14 n (shared sites)
        awk -F'\t' -v OFS='\t' -v thr=~{min_match_relatedness} -v min_n=~{min_shared_sites} '
            function expected_of(name,   s) { s = name; sub(/^(ONT|ILMN)_/, "", s); return s }
            FNR == NR { samples[++n_samples] = $1; next }
            FNR == 1 {
                print "sample_a", "sample_b", "relatedness", "n_sites" > "unexpected_matches.tsv"
                next
            }
            {
                a = $1; b = $2; r = $3 + 0; n = $14 + 0
                ea = expected_of(a); eb = expected_of(b)
                if (ea != eb && r >= thr) print a, b, r, n > "unexpected_matches.tsv"

                if (a ~ /^ILMN_/ && b ~ /^ONT_/) { t = a; a = b; b = t; t = ea; ea = eb; eb = t }
                if (a !~ /^ONT_/ || b !~ /^ILMN_/) next

                if (ea == eb) { pair_r[ea] = r; pair_n[ea] = n; ont[ea] = a; ilmn[ea] = b }
                if (!(ea in best_r) || r > best_r[ea]) { best_r[ea] = r; best_hit[ea] = b }
            }
            END {
                print "sample", "status", "pair_relatedness", "pair_n_sites", "ont_digest", "illumina_digest", "best_illumina_match", "best_relatedness"
                for (i = 1; i <= n_samples; i++) {
                    s = samples[i]
                    if (!(s in pair_r)) status = "MISSING"
                    else if (pair_r[s] >= thr && pair_n[s] >= min_n) status = "MATCH"
                    else if (pair_r[s] >= thr) status = "LOW_SITES"
                    else status = "MISMATCH"
                    print s, status, pair_r[s], pair_n[s], ont[s], ilmn[s], best_hit[s], best_r[s]
                }
            }
        ' <(printf "%s\n" ~{sep=" " sample_names}) somalier.pairs.tsv > swap_report.tsv

        cat swap_report.tsv
        cat unexpected_matches.tsv

        if awk -F'\t' 'NR > 1 && $2 != "MATCH" { bad = 1 } END { exit !bad }' swap_report.tsv
        then
            echo "FAIL" > status.txt
        else
            echo "PASS" > status.txt
        fi
    >>>

    output {
        File pairs = "working/somalier.pairs.tsv"
        File samples = "working/somalier.samples.tsv"
        File html = "working/somalier.html"
        File swap_report = "working/swap_report.tsv"
        File unexpected_matches = "working/unexpected_matches.tsv"
        String status = read_string("working/status.txt")
    }

    runtime {
        memory: memSizeGB + " GB"
        cpu: threads
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: docker_image
        # retry transient failures such as Docker Hub pull rate limits
        maxRetries: 2
    }
}
