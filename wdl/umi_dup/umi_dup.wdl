version 1.0

task UMI_dup {
    input {
        String sample_prefix
        File star_align
        Array[File] transcriptome_align = []
        Boolean emit_molecule_expression = false
        # Runtime Attributes
        Int memory
        Int disk_space
        Int ncpu
        Int preemptible
        String docker
    }

    command <<<
        set -euo pipefail

        tool_version="$(umi_tools --version)"
        if [[ "$tool_version" != "UMI-tools version: 1.1.6" ]]; then
            echo "Unexpected UMI-tools version: $tool_version" >&2
            exit 2
        fi

        echo "--- $(date "+[%b %d %H:%M:%S]") Preparing policy-eligible UMI alignments ---"
        python3 /usr/local/src/prepare_umi_bam.py \
            --input-bam "~{star_align}" \
            --output-bam eligible.bam \
            --metrics policy.tmp.json
        python3 -c 'import pysam; pysam.index("eligible.bam")'

        echo "--- $(date "+[%b %d %H:%M:%S]") Running directional UMI grouping ---"
        transcriptome_bam="~{sep="" transcriptome_align}"
        dedup_args=(
            --stdin=eligible.bam
            --stdout=/dev/stdout
            --no-sort-output
            --paired
            --extract-umi-method=tag
            --umi-tag=RX
            --method=directional
            --edit-distance-threshold=1
            --multimapping-detection-method=NH
            --unpaired-reads=discard
            --chimeric-pairs=discard
            --unmapped-reads=discard
            --random-seed=12345
            --log="~{sample_prefix}.umi_tools.log"
        )
        if [[ "~{emit_molecule_expression}" == "true" ]]; then
            if [[ "~{length(transcriptome_align)}" != "1" || ! -s "$transcriptome_bam" ]]; then
                echo "Molecule expression requires exactly one nonempty STAR transcriptome BAM" >&2
                exit 2
            fi
            umi_tools dedup "${dedup_args[@]}" | \
            python3 /usr/local/src/propagate_molecule_qnames.py \
                --representatives-bam - \
                --genomic-bam "~{star_align}" \
                --transcriptome-bam "$transcriptome_bam" \
                --genomic-output molecule.genomic.tmp.bam \
                --transcriptome-output molecule.transcriptome.tmp.bam \
                --metrics propagation.tmp.json \
                --database representatives.tmp.sqlite3 \
                --umi-length 8 \
                --representation rx_v1 \
                --container "~{docker}"
            if [[ -e representatives.tmp.sqlite3 ]]; then
                echo "Task-local representative database was not deleted" >&2
                exit 2
            fi
        else
            umi_tools dedup "${dedup_args[@]}" >/dev/null
        fi

        python3 /usr/local/src/summarize_umi_tools.py \
            --sample "~{sample_prefix}" \
            --policy policy.tmp.json \
            --log "~{sample_prefix}.umi_tools.log" \
            --metrics "~{sample_prefix}.umi_metrics.json" \
            --report "~{sample_prefix}_umi_report.txt" \
            --version "$tool_version" \
            --container "~{docker}"

        if [[ "~{emit_molecule_expression}" == "true" ]]; then
            python3 /usr/local/src/summarize_molecule_expression.py \
                --umi-metrics "~{sample_prefix}.umi_metrics.json" \
                --propagation-metrics propagation.tmp.json \
                --output molecule_expression.tmp.json
            mv molecule.genomic.tmp.bam "~{sample_prefix}.umi_molecules.genomic.bam"
            mv molecule.transcriptome.tmp.bam "~{sample_prefix}.umi_molecules.transcriptome.bam"
            mv molecule_expression.tmp.json "~{sample_prefix}.umi_molecule_expression_metrics.json"
            rm propagation.tmp.json
        fi

        rm eligible.bam eligible.bam.bai policy.tmp.json
        echo "--- $(date "+[%b %d %H:%M:%S]") Finished directional UMI grouping ---"
    >>>

    output {
        File umi_report = "${sample_prefix}_umi_report.txt"
        File umi_metrics = "${sample_prefix}.umi_metrics.json"
        File umi_log = "${sample_prefix}.umi_tools.log"
        Array[File] molecule_genomic_bam = glob("${sample_prefix}.umi_molecules.genomic.bam")
        Array[File] molecule_transcriptome_bam = glob("${sample_prefix}.umi_molecules.transcriptome.bam")
        Array[File] molecule_expression_metrics = glob("${sample_prefix}.umi_molecule_expression_metrics.json")
    }

    runtime {
        cpu: ncpu
        memory: "${memory}GB"
        disks: "local-disk ${disk_space} HDD"
        docker: docker
        preemptible: preemptible
    }

    parameter_meta {
        sample_prefix: {
            type: "id"
        }
        star_align: {
            label: "Aligned BAM File"
        }
        transcriptome_align: {
            label: "Optional STAR Transcriptome BAM File"
        }
        emit_molecule_expression: {
            label: "Emit directional UMI molecule-expression intermediates"
        }
    }

    meta {
        author: "MoTrPAC Bioinformatics Center"
        description: "Report directional UMI-tools duplicate metrics and optionally emit molecule-expression intermediates"
    }
}
