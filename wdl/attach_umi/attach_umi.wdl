version 1.0

task attachUMI {
    input {
        String SID
        File fastqr1
        File fastqr2
        File fastqi1

        # Runtime Attributes
        Int memory
        Int disk_space
        Int ncpu
        Int preemptible

        String docker
    }

    command <<<
        set -euo pipefail

        echo "--- $(date "+[%b %d %H:%M:%S]") Beginning task, making output directories ---"
        mkdir fastq_attach

        r1_tmp="fastq_attach/~{SID}_R1.fastq.gz.tmp"
        r2_tmp="fastq_attach/~{SID}_R2.fastq.gz.tmp"
        r1_pid=
        r2_pid=
        trap 'rm -f -- "$r1_tmp" "$r2_tmp"; kill $r1_pid $r2_pid 2>/dev/null || true' EXIT

        echo "--- $(date "+[%b %d %H:%M:%S]") Running attachUMI for ~{fastqr1} ---"
        (
            set -euo pipefail
            gzip -cd -- "~{fastqr1}" | gawk -v Ifq="~{fastqi1}" -f /usr/local/src/UMI_attach.awk | gzip -c > "$r1_tmp"
        ) >r1.attach.log 2>&1 &
        r1_pid=$!

        echo "--- $(date "+[%b %d %H:%M:%S]") Running attachUMI for ~{fastqr2} ---"
        (
            set -euo pipefail
            gzip -cd -- "~{fastqr2}" | gawk -v Ifq="~{fastqi1}" -f /usr/local/src/UMI_attach.awk | gzip -c > "$r2_tmp"
        ) >r2.attach.log 2>&1 &
        r2_pid=$!

        set +e
        wait "$r1_pid"
        r1_status=$?
        wait "$r2_pid"
        r2_status=$?
        set -e

        if (( r1_status != 0 || r2_status != 0 )); then
            echo "UMI attachment mate failures: R1=$r1_status R2=$r2_status" >&2
            sed 's/^/[R1] /' r1.attach.log >&2
            sed 's/^/[R2] /' r2.attach.log >&2
            exit 1
        fi

        gzip -t "$r1_tmp"
        gzip -t "$r2_tmp"

        mv -- "$r1_tmp" "fastq_attach/~{SID}_R1.fastq.gz"
        mv -- "$r2_tmp" "fastq_attach/~{SID}_R2.fastq.gz"

        trap - EXIT

        echo "--- $(date "+[%b %d %H:%M:%S]") Finished task ---"
    >>>

    output {
        File r1_umi_attached = "fastq_attach/${SID}_R1.fastq.gz"
        File r2_umi_attached = "fastq_attach/${SID}_R2.fastq.gz"
    }

    runtime {
        docker: docker
        memory: "${memory}GB"
        disks: "local-disk ${disk_space} HDD"
        cpu: ncpu
        preemptible: preemptible
    }

    parameter_meta {
        SID: {
            type: "id"
        }
        fastqr1: {
            label: "Forward End Read FASTQ File"
        }
        fastqr2: {
            label: "Reverse End Read FASTQ File"
        }
        fastqi1: {
            label: "UMI Read FASTQ File"
        }
    }

    meta {
        author: "Archana Raja"
        description: "Attach synchronized eight-base index-read UMIs to paired FASTQ headers"
    }
}
