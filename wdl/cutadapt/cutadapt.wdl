version 1.0

task Cutadapt {
    input {
        String SID
        File fastqr1
        File fastqr2

        String index_adapter
        String univ_adapter
        Int? minimumLength
        
        Int ncpu
        Int disk_space
        Int memory
        Int preemptible
        String docker
    }

    command <<<
        set -euo pipefail

        echo "--- $(date "+[%b %d %H:%M:%S]") Beginning task, making output directories ---"
        mkdir -p fastq_trim
        mkdir -p fastq_trim/tooshort

        echo "--- $(date "+[%b %d %H:%M:%S]") Running cutadapt on ~{fastqr1} and ~{fastqr2} ---"

        cutadapt \
        --cores ~{ncpu} \
        --pair-filter any \
        --compression-level 1 \
        -a ~{index_adapter} \
        -A ~{univ_adapter} \
        -o fastq_trim/~{SID}_R1.fastq.gz \
        -p fastq_trim/~{SID}_R2.fastq.gz \
        -m ~{minimumLength} \
        --too-short-output fastq_trim/tooshort/~{SID}_R1.fastq.gz \
        --too-short-paired-output fastq_trim/tooshort/~{SID}_R2.fastq.gz \
        ~{fastqr1} ~{fastqr2} > "fastq_trim/~{SID}_report.log"

        awk '
            /^Pairs written \(passing filters\):/ {
                count = $5
                gsub(",", "", count)
                if (++matches != 1 || count !~ /^[0-9]+$/) exit 2
                print count
            }
            END { if (matches != 1) exit 2 }
        ' "fastq_trim/~{SID}_report.log" > "fastq_trim/~{SID}_read_pairs.txt" || {
            echo "Could not read one passing read-pair count from the Cutadapt report" >&2
            exit 2
        }

        echo "--- $(date "+[%b %d %H:%M:%S]") Cutadapt done, task complete ---"
    >>>

    output {
        File fastq_trimmed_R1 = "fastq_trim/${SID}_R1.fastq.gz"
        File fastq_trimmed_R2 = "fastq_trim/${SID}_R2.fastq.gz"
        File report = "fastq_trim/${SID}_report.log"
        Int read_pairs = read_int("fastq_trim/${SID}_read_pairs.txt")
        File tooShortOutput = "fastq_trim/tooshort/${SID}_R1.fastq.gz"
        File tooShortPairedOutput = "fastq_trim/tooshort/${SID}_R2.fastq.gz"
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
    }

    runtime {
        cpu: ncpu
        memory: "${memory}GB"
        disks: "local-disk ${disk_space} HDD"
        docker: docker
        preemptible: preemptible
    }
}
