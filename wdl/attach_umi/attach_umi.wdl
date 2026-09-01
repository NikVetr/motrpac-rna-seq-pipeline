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
        trap 'rm -f -- "$r1_tmp" "$r2_tmp"' EXIT

        echo "--- $(date "+[%b %d %H:%M:%S]") Running attachUMI for ~{fastqr1} ---"
        gzip -cd -- "~{fastqr1}" | gawk -v Ifq="~{fastqi1}" -f /usr/local/src/UMI_attach.awk | gzip -c > "$r1_tmp"
        gzip -t "$r1_tmp"

        echo "--- $(date "+[%b %d %H:%M:%S]") Running attachUMI for ~{fastqr2} ---"
        gzip -cd -- "~{fastqr2}" | gawk -v Ifq="~{fastqi1}" -f /usr/local/src/UMI_attach.awk | gzip -c > "$r2_tmp"
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
