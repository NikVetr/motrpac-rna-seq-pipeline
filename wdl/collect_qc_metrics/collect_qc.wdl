version 1.0

task rnaseqQC {
    input {
        String SID
        File? fastqc_pretrim
        File? fastqc_posttrim
        String pretrim_r1_filename
        String pretrim_r2_filename
        String posttrim_r1_filename
        String posttrim_r2_filename
        File cutadapt_report
        File? mapped_report
        File? rRNA_report
        File? globin_report
        File? phix_report
        File star_log
        File? markduplicates_metrics
        File? rnaseq_metrics
        File? umi_report

        Int memory
        Int disk_space
        Int ncpu
        Int preemptible
        String docker
    }

    command <<<
        set -euo pipefail
        echo "--- $(date "+[%b %d %H:%M:%S]") Collecting native RNA-seq QC reports ---"

        python3 /usr/local/src/rnaseq_qc.py \
            --sample "~{SID}" \
            ~{"--fastqc-pretrim \"" + fastqc_pretrim + "\""} \
            ~{"--fastqc-posttrim \"" + fastqc_posttrim + "\""} \
            --pretrim-r1-filename "~{pretrim_r1_filename}" \
            --pretrim-r2-filename "~{pretrim_r2_filename}" \
            --posttrim-r1-filename "~{posttrim_r1_filename}" \
            --posttrim-r2-filename "~{posttrim_r2_filename}" \
            --cutadapt-report "~{cutadapt_report}" \
            ~{"--mapped-report \"" + mapped_report + "\""} \
            ~{"--rrna-report \"" + rRNA_report + "\""} \
            ~{"--globin-report \"" + globin_report + "\""} \
            ~{"--phix-report \"" + phix_report + "\""} \
            --star-log "~{star_log}" \
            ~{"--markduplicates-metrics \"" + markduplicates_metrics + "\""} \
            ~{"--rnaseq-metrics \"" + rnaseq_metrics + "\""} \
            ~{"--umi-report \"" + umi_report + "\""} \
            --output "~{SID}_qc_info.csv"

        echo "--- $(date "+[%b %d %H:%M:%S]") Finished native RNA-seq QC collection ---"
    >>>

    output {
        File rnaseq_report = "${SID}_qc_info.csv"
    }

    runtime {
        cpu: ncpu
        memory: "${memory}GB"
        disks: "local-disk ${disk_space} HDD"
        docker: docker
        preemptible: preemptible
    }

    parameter_meta {
        SID: {
            type: "id"
        }
        fastqc_pretrim: {
            label: "Pre-trim FastQC archive"
        }
        fastqc_posttrim: {
            label: "Post-trim FastQC archive"
        }
        umi_report: {
            label: "Optional UMI duplicate-rate report"
        }
    }

    meta {
        author: "Archana Raja"
        description: "Collect published QC fields directly from native task reports"
    }
}
