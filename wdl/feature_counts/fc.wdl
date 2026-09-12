version 1.0

task feature_counts {
    input {
        String SID
        File input_bam
        File gtf_file

        Int memory
        Int disk_space
        Int ncpu
        Boolean use_e2 = false
        Int preemptible
        String docker
    }

    Int e2_cpu = 2 * ceil(if ncpu / 2.0 > memory / 16.0 then ncpu / 2.0 else memory / 16.0)

    command <<<
        set -euo pipefail
        echo "--- $(date "+[%b %d %H:%M:%S]") Beginning task, running featurecounts ---"
        featureCounts \
            -T ~{ncpu} \
            -a ~{gtf_file} \
            -o ~{SID}.out \
            -p \
            --countReadPairs \
            -s 1 \
            -M \
            --fraction \
            ~{input_bam}

        echo "$(date "+[%b %d %H:%M:%S]") Finished featurecounts"
        ls -ltr

        echo "--- $(date "+[%b %d %H:%M:%S]") Finished task ---"
    >>>

    output {
        File fc_out = "${SID}.out"
        File fc_summary = "${SID}.out.summary"
    }

    runtime {
        docker: "${docker}"
        memory: "${memory}GB"
        disks: "local-disk ${disk_space} HDD"
        cpu: "${ncpu}"
        gcp: if use_e2 then object { predefinedMachineType: "e2-custom-${e2_cpu}-${memory * 1024}" } else object {}
        preemptible: preemptible
    }

    parameter_meta {
        SID: {
            type: "id"
        }
        input_bam: {
            label: "Input BAM File"
        }
        gtf_file: {
            label: "GTF-Format Annotation File"
        }
    }

    meta {
        author: "Archana Raja"
        description: "Forward-stranded paired-fragment gene counting"
    }
}
