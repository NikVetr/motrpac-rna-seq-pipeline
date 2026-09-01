version 1.0

task fastQC {
    input {
        File fastqr1
        File fastqr2

        String outdir

        Int memory
        Int disk_space
        Int ncpu
        Int preemptible
        String docker
    }

    command <<<
        set -euo pipefail

        echo "--- $(date "+[%b %d %H:%M:%S]") Beginning task, running mate-level FastQC ---"

        if (( ~{ncpu} < 2 )); then
            echo "parallel FastQC requires at least two CPUs; received ~{ncpu}" >&2
            exit 2
        fi

        mkdir -p -- "~{outdir}" r1_fastqc r2_fastqc

        fastqc -o r1_fastqc -- "~{fastqr1}" >r1.fastqc.log 2>&1 &
        r1_pid=$!
        fastqc -o r2_fastqc -- "~{fastqr2}" >r2.fastqc.log 2>&1 &
        r2_pid=$!

        set +e
        wait "$r1_pid"
        r1_status=$?
        wait "$r2_pid"
        r2_status=$?
        set -e

        if (( r1_status != 0 || r2_status != 0 )); then
            echo "FastQC mate failures: R1=$r1_status R2=$r2_status" >&2
            sed 's/^/[R1] /' r1.fastqc.log >&2
            sed 's/^/[R2] /' r2.fastqc.log >&2
            exit 1
        fi

        shopt -s nullglob
        r1_zip=(r1_fastqc/*_fastqc.zip)
        r1_html=(r1_fastqc/*_fastqc.html)
        r2_zip=(r2_fastqc/*_fastqc.zip)
        r2_html=(r2_fastqc/*_fastqc.html)
        if (( ${#r1_zip[@]} != 1 || ${#r1_html[@]} != 1 ||
              ${#r2_zip[@]} != 1 || ${#r2_html[@]} != 1 )); then
            echo "FastQC did not emit exactly one ZIP and HTML report per mate" >&2
            exit 1
        fi
        if [[ "${r1_zip[0]##*/}" == "${r2_zip[0]##*/}" ||
              "${r1_html[0]##*/}" == "${r2_html[0]##*/}" ]]; then
            echo "FastQC report names collide between mates" >&2
            exit 1
        fi

        mv -- "${r1_zip[0]}" "${r1_html[0]}" "${r2_zip[0]}" "${r2_html[0]}" "~{outdir}/"
        tar -czf "~{outdir}.tar.gz" "./~{outdir}"
        echo "--- $(date "+[%b %d %H:%M:%S]") Finished mate-level FastQC, task complete ---"
    >>>

    output {
        File fastQC_report = "${outdir}.tar.gz"
    }

    parameter_meta {
        fastqr1: {
            label: "Forward End Read FASTQ File"
        }
        fastqr2: {
            label: "Reverse End Read FASTQ File"
        }
    }

    runtime {
        cpu: ncpu
        memory: "${memory}GB"
        disks: "local-disk ${disk_space} HDD"
        docker: docker
        preemptible: preemptible
    }

    meta {
        description: "Run mate-level FastQC processes concurrently with explicit status and report validation"
    }
}
