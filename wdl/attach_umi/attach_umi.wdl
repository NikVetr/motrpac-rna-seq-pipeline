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
        Boolean use_e2 = false
        Int preemptible

        String docker
    }

    Int e2_cpu = 2 * ceil(if ncpu / 2.0 > memory / 16.0 then ncpu / 2.0 else memory / 16.0)

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
        gcp: if use_e2 then object { predefinedMachineType: "e2-custom-${e2_cpu}-${memory * 1024}" } else object {}
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

task trimIndexRead {
    input {
        String SID
        File fastqi1
        Boolean use_e2 = false
        Int preemptible
        String docker
    }

    Int disk_gb = ceil(10.0 + 3.0 * size(fastqi1, "GiB"))

    command <<<
        set -euo pipefail
        # Eight-base I1 passes through untouched; attachUMI validates every record.
        # Every nine-base record is trimmed; validate the extra-A layout over the whole file.
        gawk -v src="~{fastqi1}" -v out="~{SID}_I1.trimmed.fastq.gz" '
            function fail(message) { print "trimIndexRead: " message > "/dev/stderr"; exit 1 }
            function read_record() {
                if ((reader | getline header) <= 0) return 0
                if ((reader | getline seq) <= 0 || (reader | getline plus) <= 0 || (reader | getline qual) <= 0 ||
                    substr(header, 1, 1) != "@" || substr(plus, 1, 1) != "+") fail("malformed record " n + 1)
                return 1
            }
            BEGIN {
                reader = "gzip -cd -- " src
                if (!read_record()) fail("index FASTQ is empty")
                if (length(seq) == 8) exit
                if (length(seq) != 9) fail("unexpected UMI length " length(seq))
                writer = "gzip -c > " out
                do {
                    if (length(seq) != 9 || length(qual) != 9) fail("record " n + 1 " is not nine bases")
                    print header "\n" substr(seq, 1, 8) "\n+\n" substr(qual, 1, 8) | writer
                    if (substr(seq, 9, 1) == "A") adenines++
                    n++
                } while (read_record())
                if (close(reader) || close(writer)) fail("index decompression or compression failed")
                printf "trimIndexRead: records=%d ninth_A_count=%d ninth_A_fraction=%.6f\n", n, adenines, adenines / n
                if (adenines < 0.9 * n) fail("ninth-base A fraction is below 90%; not the expected extra-A layout")
            }'
    >>>

    output {
        Array[File] trimmed_index = glob("${SID}_I1.trimmed.fastq.gz")
    }

    runtime {
        docker: docker
        cpu: 1
        memory: "2GB"
        gcp: if use_e2 then object { predefinedMachineType: "e2-custom-2-2048" } else object {}
        disks: "local-disk ${disk_gb} HDD"
        preemptible: preemptible
    }

    meta {
        description: "Optionally trim a fixed ninth base from nine-base I1 reads carrying an eight-base UMI"
    }
}
