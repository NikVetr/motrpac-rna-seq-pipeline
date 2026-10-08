version 1.0

import "merge_results/merge_results.wdl" as final_merge
import "merge_results/merge_isoforms.wdl" as isoform_merge

# Merge per-sample results saved by one or more rnaseq_pipeline runs, without recomputation.
# Inputs are normally written by scripts/make_merge_inputs.py.
workflow rnaseq_merge {
    input {
        # Completed samples, in matrix column order.
        Array[String]+ sample_prefix
        # Header plus one row per cohort sample, completed or failed, in cohort order.
        Array[Array[String]] expression_metadata_rows
        Array[String] failed_samples = []
        Array[File] rsem_gene_results
        Array[File] rsem_isoform_results
        Array[File] feature_counts_files
        Array[File] qc_report_files
        Array[File] qc_diagnostics
        Array[Array[String]] sample_size_rows = []
        String output_report_name
        Int merge_results_ncpu = 1
        Int merge_results_ramGB = 4
        Int merge_results_disk = 10
        Int num_preemptible_attempts = 0
        String merge_results_docker
    }

    call final_merge.merge_results {
        input:
            sample_prefix=sample_prefix,
            output_report_name=output_report_name,
            rsem_files=rsem_gene_results,
            feature_counts_files=feature_counts_files,
            qc_report_files=qc_report_files,
            expression_metadata_rows=expression_metadata_rows,
            qc_diagnostics=qc_diagnostics,
            sample_size_rows=sample_size_rows,
            failed_samples=failed_samples,
            ncpu=merge_results_ncpu,
            memory=merge_results_ramGB,
            disk_space=merge_results_disk,
            preemptible=num_preemptible_attempts,
            docker=merge_results_docker
    }

    call isoform_merge.merge_isoforms {
        input:
            sample_prefix=sample_prefix,
            rsem_files=rsem_isoform_results,
            rsem_gene_files=rsem_gene_results,
            ncpu=merge_results_ncpu,
            memory=merge_results_ramGB,
            disk_space=merge_results_disk,
            preemptible=num_preemptible_attempts,
            docker=merge_results_docker
    }

    output {
        File rsem_genes_count = merge_results.rsem_genes_count
        File rsem_genes_tpm = merge_results.rsem_genes_tpm
        File rsem_genes_fpkm = merge_results.rsem_genes_fpkm
        File rsem_genes_effective_length = merge_results.rsem_genes_effective_length
        File rsem_isoforms_count = merge_isoforms.rsem_isoforms_count
        File rsem_isoforms_tpm = merge_isoforms.rsem_isoforms_tpm
        File rsem_isoforms_fpkm = merge_isoforms.rsem_isoforms_fpkm
        File rsem_isoforms_effective_length = merge_isoforms.rsem_isoforms_effective_length
        File rsem_transcripts = merge_isoforms.rsem_transcripts
        File feature_counts_file = merge_results.feature_counts
        File qc_report_file = merge_results.qc_report
        File expression_metadata = merge_results.expression_metadata
    }
}
