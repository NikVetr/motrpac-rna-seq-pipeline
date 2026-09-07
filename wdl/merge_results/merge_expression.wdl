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
    # Four GiB per 75 libraries gives the validated 16-GiB request at 297.
    Int inferred_merge_memory_gb = 4 * ceil(length(sample_prefix) / 75.0)
    Int effective_merge_memory_gb =
        if memory > inferred_merge_memory_gb then memory else inferred_merge_memory_gb
    Float merge_input_gib = size(rsem_files, "GiB") + size(feature_counts_files, "GiB")
    # Retain the cohort-validated allowance for localization, outputs, and scratch.
    Int inferred_merge_scratch_gb = ceil(3.0 * merge_input_gib + 10.0)
    Int effective_merge_scratch_gb =
        if disk_space > inferred_merge_scratch_gb then disk_space else inferred_merge_scratch_gb

    command <<<
        set -euo pipefail
        mkdir rsem_files feature_counts_files
        while IFS= read -r input; do
            ln -s -- "$(realpath -- "$input")" rsem_files/
        done < "~{write_lines(rsem_files)}"
        while IFS= read -r input; do
            ln -s -- "$(realpath -- "$input")" feature_counts_files/
        done < "~{write_lines(feature_counts_files)}"

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
        memory: "${effective_merge_memory_gb}GB"
        disks: "local-disk ${effective_merge_scratch_gb} HDD"
        docker: docker
        preemptible: preemptible
    }

    meta {
        author: "MoTrPAC Bioinformatics Center"
        description: "Merge optional secondary expression results"
    }
}
