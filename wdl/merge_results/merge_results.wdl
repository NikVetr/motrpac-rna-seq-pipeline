version 1.0

task merge_results {
    input {
        Array[String]+ sample_prefix
        Array[File] rsem_files
        Array[File] feature_counts_files
        Array[File] qc_report_files
        String output_report_name

        Int memory
        Int disk_space
        Int ncpu
        Int preemptible
        String docker
    }

    File sample_order = write_lines(sample_prefix)
    Float merge_input_gib = size(rsem_files, "GiB") + size(feature_counts_files, "GiB") + size(qc_report_files, "GiB")
    # Allow for localized inputs, task-local copies, merged outputs, and scratch.
    Int inferred_merge_scratch_gb = ceil(3.0 * merge_input_gib + 10.0)
    Int effective_merge_scratch_gb =
        if disk_space > inferred_merge_scratch_gb then disk_space else inferred_merge_scratch_gb

    command <<<
        set -eou pipefail
        echo "--- $(date "+[%b %d %H:%M:%S]") Beginning task, copying files ---"

        mkdir -p rsem_files
        mkdir -p qc_report_files
        mkdir -p feature_counts_files

        cp ~{sep=" " rsem_files} rsem_files/
        cp ~{sep=" " feature_counts_files} feature_counts_files/
        cp ~{sep=" " qc_report_files} qc_report_files/

        echo "--- $(date "+[%b %d %H:%M:%S]") Merging RSEM results ---"
        python3 /usr/local/src/merge_rsem.py \
            --rsem-dir rsem_files \
            --sample-order ~{sample_order}

        echo "--- $(date "+[%b %d %H:%M:%S]") Finished merging RSEM results, consolidating QC reports ---"
        python3 /usr/local/src/consolidate_qc_report.py \
            --qc-dir qc_report_files \
            --sample-order ~{sample_order} \
            --output-name ~{output_report_name}.csv

        echo "--- $(date "+[%b %d %H:%M:%S]") Finished merging consolidating QC reports, merging feature counts ---"
        python3 /usr/local/src/merge_fc.py \
            --fc-dir feature_counts_files \
            --sample-order ~{sample_order}

        echo "--- $(date "+[%b %d %H:%M:%S]") Finished merging feature counts, finished task  ---"
    >>>

    output {
        File rsem_genes_count = "rsem_genes_count.txt"
        File rsem_genes_tpm = "rsem_genes_tpm.txt"
        File rsem_genes_fpkm = "rsem_genes_fpkm.txt"
        File feature_counts = "featureCounts.txt"
        File qc_report = "${output_report_name}.csv"
    }

    parameter_meta {
        sample_prefix: {
            type: "id"
        }
        rsem_files: {
            label: "RSEM Gene Count Files"
        }
        feature_counts_files: {
            label: "Feature Counts Files"
        }
        qc_report_files: {
            label: "QC Report Files"
        }
    }

    runtime {
        cpu: ncpu
        memory: "${memory}GB"
        disks: "local-disk ${effective_merge_scratch_gb} HDD"
        docker: docker
        preemptible: preemptible
    }
}
