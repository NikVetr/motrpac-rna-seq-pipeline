version 1.0

task rsem {
    input {
        String SID
        File transcriptome_bam
        File rsem_reference
        String reference_release = "unspecified"
        Boolean umi_deduplicated = false

        Int memory
        Int disk_space
        Int ncpu
        Boolean use_e2 = false
        Int preemptible
        String docker
    }

    # Release-specific working-memory buffers; size the BAM entering this call.
    Float input_gib = size(transcriptome_bam, "GiB")
    Int inferred_memory = if reference_release == "rn8_v116" && umi_deduplicated then ceil(2.0 + 4.0 * input_gib)
                          else if reference_release == "rn8_v116" then ceil(if input_gib <= 4.375 then 18.0 else 0.5 + 4.0 * input_gib)
                          else if reference_release == "gencode_v50" && umi_deduplicated then 2 * ceil((14.0 + 1.6 * input_gib) / 2.0)
                          else 4 * ceil((16.0 + 2.0 * input_gib) / 4.0)
    Int inferred_scratch_gb = ceil(10.0 + 4.0 * input_gib)
    Int effective_memory = if memory > inferred_memory then memory else inferred_memory
    Int effective_scratch_gb = if disk_space > inferred_scratch_gb then disk_space else inferred_scratch_gb
    Int e2_cpu = 2 * ceil(if ncpu / 2.0 > effective_memory / 16.0 then ncpu / 2.0 else effective_memory / 16.0)

    command <<<
        set -euo pipefail
        mkdir rsem_reference
        echo "$(date "+[%b %d %H:%M:%S]") Extracting rsem_reference"
        tar -xzvf ~{rsem_reference} -C rsem_reference --strip-components=1
        echo "$(date "+[%b %d %H:%M:%S]") Done tar"

        cd rsem_reference
        echo "--- Running: ls --- "
        ls
        echo "--- $(date "+[%b %d %H:%M:%S]") Running: rsem-calculate-expression --- "
        rsem-calculate-expression \
            -p ~{ncpu} \
            --bam \
            --paired-end \
            --no-bam-output \
            --forward-prob 1 \
            --seed 12345 \
            ~{transcriptome_bam} \
            rsem_reference \
            ~{SID} 2>&1 | tee ~{SID}.rsem.log
        echo "--- $(date "+[%b %d %H:%M:%S]") Done: rsem-calculate-expression --- "
        ls
        echo "--- $(date "+[%b %d %H:%M:%S]") Finished task --- "
    >>>

    output {
        File genes = "rsem_reference/${SID}.genes.results"
        File isoforms = "rsem_reference/${SID}.isoforms.results"
        File stat_cnt = "rsem_reference/${SID}.stat/${SID}.cnt"
        File stat_model = "rsem_reference/${SID}.stat/${SID}.model"
        File stat_theta = "rsem_reference/${SID}.stat/${SID}.theta"
        File log = "rsem_reference/${SID}.rsem.log"
        File convergence = "rsem_reference/${SID}.rsem_convergence.tsv"
        File gene_convergence = "rsem_reference/${SID}.rsem_gene_convergence.tsv"
    }

    runtime {
        cpu: ncpu
        gcp: if use_e2 then object { predefinedMachineType: "e2-custom-${e2_cpu}-${effective_memory * 1024}" } else object {}
        memory: "${effective_memory}GB"
        disks: "local-disk ${effective_scratch_gb} HDD"
        docker: docker
        preemptible: preemptible
    }

    parameter_meta {
        SID: {
            type: "id"
        }
        transcriptome_bam: {
            label: "Aligned Transcriptome BAM File"
        }
        rsem_reference: {
            label: "RSEM Genome Reference File"
        }
    }

    meta {
        author: "Archana Raja"
        description: "Forward-stranded paired-end transcript quantification"
    }
}
