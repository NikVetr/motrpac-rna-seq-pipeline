version 1.0

task merge_isoforms {
    input {
        Array[String]+ sample_prefix
        Array[File] rsem_files
        # Gene results from the same calls; each gene must equal the sum of its isoforms.
        Array[File] rsem_gene_files
        Int memory
        Int disk_space
        Int ncpu
        Int preemptible
        String docker
    }

    File sample_order = write_lines(sample_prefix)
    Int inferred_scratch_gb = ceil(3.0 * (size(rsem_files, "GiB") + size(rsem_gene_files, "GiB")) + 10.0)
    Int effective_scratch_gb = if disk_space > inferred_scratch_gb then disk_space else inferred_scratch_gb
    # Streaming rows avoids holding the transcript-by-sample matrix in memory.
    Int effective_memory_gb = if memory > 4 then memory else 4

    command <<<
        set -euo pipefail
        mkdir rsem_files rsem_gene_files
        while IFS= read -r input; do
            ln -s -- "$(realpath -- "$input")" rsem_files/
        done < "~{write_lines(rsem_files)}"
        while IFS= read -r input; do
            ln -s -- "$(realpath -- "$input")" rsem_gene_files/
        done < "~{write_lines(rsem_gene_files)}"
        python3 /usr/local/src/merge_rsem.py \
            --rsem-dir rsem_files --sample-order "~{sample_order}" --feature-level isoforms \
            --gene-dir rsem_gene_files
    >>>

    output {
        File rsem_isoforms_count = "rsem_isoforms_count.txt"
        File rsem_isoforms_tpm = "rsem_isoforms_tpm.txt"
        File rsem_isoforms_fpkm = "rsem_isoforms_fpkm.txt"
        File rsem_isoforms_effective_length = "rsem_isoforms_effective_length.txt"
        File rsem_transcripts = "rsem_transcripts.tsv"
    }

    runtime {
        cpu: ncpu
        memory: "${effective_memory_gb}GB"
        disks: "local-disk ${effective_scratch_gb} HDD"
        docker: docker
        preemptible: preemptible
    }
}
