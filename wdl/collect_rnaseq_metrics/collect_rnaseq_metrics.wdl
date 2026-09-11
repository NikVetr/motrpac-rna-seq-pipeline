version 1.0

#Change MINIMUM_LENGTH=50 RRNA_FRAGMENT_PERCENTAGE=0.3, present in the shell script was missing in the MOP

task collectrnaseqmetrics {
    input {
        String SID
        File input_bam
        File ref_flat
        
        Int memory
        Int disk_space
        Int ncpu
        Boolean prefer_predefined_n1 = false
        Int preemptible
        String docker
    }

    command <<<
        set -euo pipefail
        echo "--- $(date "+[%b %d %H:%M:%S]") Beginning task, making output directories ---"
        mkdir -p qc53
        mkdir -p qc53/log

        echo "--- $(date "+[%b %d %H:%M:%S]") Running Picard collect metrics ---"
        picard -Xmx~{memory}g CollectRnaSeqMetrics \
            I=~{input_bam} \
            O=qc53/~{SID}.RNA_Metrics \
            REF_FLAT=~{ref_flat} \
            STRAND=FIRST_READ_TRANSCRIPTION_STRAND \
            MINIMUM_LENGTH=50 \
            RRNA_FRAGMENT_PERCENTAGE=0.3 >& qc53/log/~{SID}.log

        ls -la qc53

        echo "--- $(date "+[%b %d %H:%M:%S]") Task complete ---"
    >>>

    output {
        File rnaseqmetrics = "qc53/${SID}.RNA_Metrics"
        File log = "qc53/log/${SID}.log"
    }

    runtime {
        cpu: ncpu
        # Same-family upgrade, cheaper in both us-west2 markets; other sizes stay custom.
        gcp: if prefer_predefined_n1 && ncpu == 2 && memory == 12
            then object { predefinedMachineType: "n1-highmem-2" } else object {}
        memory: "${memory}GB"
        disks: "local-disk ${disk_space} HDD"
        docker: docker
        preemptible: preemptible
    }

    parameter_meta {
        SID: {
            type: "id"
        }
        input_bam: {
            label: "Aligned BAM file"
        }
        ref_flat: {
            label: "RNA Transcript Reference File"
        }
    }

    meta {
        author: "Archana Raja"
    }
}
