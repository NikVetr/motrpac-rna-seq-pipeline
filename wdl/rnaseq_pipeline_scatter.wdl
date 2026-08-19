version 1.0

import "umi_dup/umi_dup.wdl" as umi_dup

workflow rnaseq_pipeline {

    meta {
        task_labels: {
            udup: {
                task_name: "UMI Duplication",
                description: "Runs nudup.py python2 script to get an estimation of the PCR duplicates rate from STAR-aligned BAM file"
            }
        }
    }

    parameter_meta {
        sample_prefix: {
            type: "id"
        }
    }

    input {
        # Inputs
        Array[File]+ star_align_bam
        Array[String]+ sample_prefix

        # Runtime controls
        Int num_preemptible_attempts = 0

        # UMI Duplication Parameters
        Int umi_dup_ncpu
        Int umi_dup_ramGB
        Int umi_dup_disk
        String umi_dup_docker
    }

    scatter (i in range(length(star_align_bam))) {
        call umi_dup.UMI_dup as udup {
            input:
            # Inputs
                sample_prefix=sample_prefix[i],
                star_align=star_align_bam[i],
            # Runtime Parameters
                ncpu=umi_dup_ncpu,
                memory=umi_dup_ramGB,
                disk_space=umi_dup_disk,
                preemptible=num_preemptible_attempts,
                docker=umi_dup_docker
        }
    }

    output {
        Array[File] umi_dup_logs = udup.umi_dup_out
        Array[File] umi_outputs = udup.umi_out
        Array[File] umi_reports = udup.umi_report
    }
}
