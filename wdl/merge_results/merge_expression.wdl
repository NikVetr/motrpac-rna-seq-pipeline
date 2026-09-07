version 1.0

task merge_expression {
    input {
        Array[String]+ sample_prefix
        String output_prefix
        Array[File] rsem_files
        Array[File] feature_counts_files

        Int memory
        Int disk_space
        Int ncpu
        Int preemptible
        String docker
    }

    File sample_order = write_lines(sample_prefix)
    Float merge_input_gib = size(rsem_files, "GiB") + size(feature_counts_files, "GiB")
    # Allow for localized inputs, task-local copies, merged outputs, and scratch.
    Int inferred_merge_scratch_gb = ceil(3.0 * merge_input_gib + 10.0)
    Int effective_merge_scratch_gb =
        if disk_space > inferred_merge_scratch_gb then disk_space else inferred_merge_scratch_gb

    command <<<
        set -euo pipefail
        mkdir rsem_files feature_counts_files
        cp ~{sep=" " rsem_files} rsem_files/
        cp ~{sep=" " feature_counts_files} feature_counts_files/

        python3 /usr/local/src/merge_rsem.py \
            --rsem-dir rsem_files \
            --sample-order ~{sample_order}
        python3 /usr/local/src/merge_fc.py \
            --fc-dir feature_counts_files \
            --sample-order ~{sample_order}

        mv rsem_genes_count.txt "~{output_prefix}_rsem_genes_count.txt"
        mv rsem_genes_tpm.txt "~{output_prefix}_rsem_genes_tpm.txt"
        mv rsem_genes_fpkm.txt "~{output_prefix}_rsem_genes_fpkm.txt"
        mv featureCounts.txt "~{output_prefix}_featureCounts.txt"
    >>>

    output {
        File rsem_genes_count = "${output_prefix}_rsem_genes_count.txt"
        File rsem_genes_tpm = "${output_prefix}_rsem_genes_tpm.txt"
        File rsem_genes_fpkm = "${output_prefix}_rsem_genes_fpkm.txt"
        File feature_counts = "${output_prefix}_featureCounts.txt"
    }

    runtime {
        cpu: ncpu
        memory: "${memory}GB"
        disks: "local-disk ${effective_merge_scratch_gb} HDD"
        docker: docker
        preemptible: preemptible
    }

    meta {
        author: "MoTrPAC Bioinformatics Center"
        description: "Merge optional secondary expression results"
    }
}
