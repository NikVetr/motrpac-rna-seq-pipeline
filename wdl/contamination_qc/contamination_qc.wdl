version 1.0

task contamination_qc {
    input {
        String SID
        File fastqr1
        File fastqr2
        File cutadapt_report
        File globin_genome_dir_tar
        File rrna_genome_dir_tar
        File phix_genome_dir_tar
        Int sample_pairs

        Int memory
        Int disk_space
        Int ncpu
        Int preemptible
        String docker
    }

    command <<<
        set -euo pipefail

        python3 - \
            "~{fastqr1}" \
            "~{fastqr2}" \
            "~{cutadapt_report}" \
            "~{sample_pairs}" \
            "motrpac-contamination-qc-v1" \
            "~{SID}.sampled_R1.fastq.gz" \
            "~{SID}.sampled_R2.fastq.gz" \
            "~{SID}_contamination_sampling.json" <<'PYTHON'
import gzip
import hashlib
import json
import re
import sys
from pathlib import Path


def read_record(handle, mate, record_number):
    lines = [handle.readline() for _ in range(4)]
    if not any(lines):
        return None
    if not all(lines):
        raise ValueError(
            "{} FASTQ is truncated at record {}".format(mate, record_number)
        )
    header, sequence, plus, quality = lines
    if not header.startswith(b"@") or not plus.startswith(b"+"):
        raise ValueError(
            "{} FASTQ has a malformed record {}".format(mate, record_number)
        )
    sequence = sequence.rstrip(b"\r\n")
    quality = quality.rstrip(b"\r\n")
    if not sequence or len(sequence) != len(quality):
        raise ValueError(
            "{} FASTQ has unequal sequence/quality lengths at record {}".format(
                mate, record_number
            )
        )
    header_fields = header[1:].split()
    token = header_fields[0] if header_fields else b""
    if not token:
        raise ValueError(
            "{} FASTQ has an empty read name at record {}".format(
                mate, record_number
            )
        )
    if token.endswith((b"/1", b"/2")):
        expected_suffix = b"/" + mate.encode("ascii")
        if not token.endswith(expected_suffix):
            raise ValueError(
                "{} FASTQ has the wrong mate suffix at record {}".format(
                    mate, record_number
                )
            )
        token = token[:-2]
    if len(header_fields) > 1:
        read_number = header_fields[1].split(b":", 1)[0]
        if read_number in (b"1", b"2") and read_number != mate.encode("ascii"):
            raise ValueError(
                "{} FASTQ has the wrong read-number field at record {}".format(
                    mate, record_number
                )
            )
    return lines, token


def paired_records(r1_path, r2_path):
    with gzip.open(r1_path, "rb") as r1_handle, gzip.open(r2_path, "rb") as r2_handle:
        record_number = 1
        while True:
            r1_record = read_record(r1_handle, "1", record_number)
            r2_record = read_record(r2_handle, "2", record_number)
            if r1_record is None and r2_record is None:
                return
            if r1_record is None or r2_record is None:
                raise ValueError(
                    "paired FASTQs contain different record counts at record {}".format(
                        record_number
                    )
                )
            r1_lines, r1_name = r1_record
            r2_lines, r2_name = r2_record
            if r1_name != r2_name:
                raise ValueError(
                    "paired FASTQs are unsynchronized at record {}: {!r} != {!r}".format(
                        record_number, r1_name, r2_name
                    )
                )
            yield r1_lines, r2_lines, r1_name
            record_number += 1


def cutadapt_pairs(path):
    text = Path(path).read_text(encoding="utf-8")
    matches = re.findall(
        r"^Pairs written \(passing filters\):\s+([0-9][0-9,]*)\s+\([^)]+\)\s*$",
        text,
        flags=re.MULTILINE,
    )
    if len(matches) != 1:
        raise ValueError(
            "Cutadapt report must contain one Pairs written (passing filters) value"
        )
    return int(matches[0].replace(",", ""))


def select_indices(population_size, sample_pairs, seed):
    selected = set()
    seed_bytes = seed.encode("utf-8")
    counter = 0

    def randbelow(bound):
        nonlocal counter
        range_size = 1 << 256
        limit = range_size - (range_size % bound)
        while True:
            payload = seed_bytes + b"\0" + counter.to_bytes(16, "big")
            counter += 1
            value = int.from_bytes(hashlib.sha256(payload).digest(), "big")
            if value < limit:
                return value % bound

    selected_count = min(sample_pairs, population_size)
    for upper in range(population_size - selected_count, population_size):
        candidate = randbelow(upper + 1)
        selected.add(upper if candidate in selected else candidate)
    if len(selected) != selected_count:
        raise AssertionError("deterministic ordinal sampler selected the wrong count")
    return selected


def write_sample(r1_path, r2_path, selected, expected_pairs, r1_output, r2_output):
    selected_names = hashlib.sha256()
    written = 0
    with open(r1_output, "wb") as r1_raw, open(r2_output, "wb") as r2_raw:
        with gzip.GzipFile(
            filename="", mode="wb", compresslevel=1, mtime=0, fileobj=r1_raw
        ) as r1_gzip, gzip.GzipFile(
            filename="", mode="wb", compresslevel=1, mtime=0, fileobj=r2_raw
        ) as r2_gzip:
            observed_pairs = 0
            for index, (r1_lines, r2_lines, read_name) in enumerate(
                paired_records(r1_path, r2_path)
            ):
                observed_pairs += 1
                if index not in selected:
                    continue
                r1_gzip.writelines(r1_lines)
                r2_gzip.writelines(r2_lines)
                selected_names.update(read_name)
                selected_names.update(b"\n")
                written += 1
    if observed_pairs != expected_pairs:
        raise ValueError(
            "Cutadapt reports {} pairs, but the FASTQs contain {}".format(
                expected_pairs, observed_pairs
            )
        )
    if written != len(selected):
        raise ValueError(
            "selected {} records but wrote {}".format(len(selected), written)
        )
    return selected_names.hexdigest()


def main():
    if len(sys.argv) != 9:
        raise ValueError(
            "expected R1, R2, Cutadapt report, N, seed, output R1, output R2, and manifest"
        )
    r1_path, r2_path, cutadapt_path = sys.argv[1:4]
    try:
        sample_pairs = int(sys.argv[4])
    except ValueError as error:
        raise ValueError("sample-pair count must be an integer") from error
    if sample_pairs < 0:
        raise ValueError("sample-pair count must be nonnegative")
    seed, r1_output, r2_output, manifest_path = sys.argv[5:]
    if sample_pairs == 0:
        input_pairs = cutadapt_pairs(cutadapt_path)
        manifest = {
            "algorithm": "none-full-depth",
            "cutadapt_report": Path(cutadapt_path).name,
            "input_fastq_1": Path(r1_path).name,
            "input_fastq_2": Path(r2_path).name,
            "input_pairs": input_pairs,
            "requested_pairs": 0,
            "seed": None,
            "selected_name_sha256": None,
            "selected_pairs": input_pairs,
            "used_full_input": True,
        }
        Path(manifest_path).write_text(
            json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )
        return
    input_pairs = cutadapt_pairs(cutadapt_path)
    selected = select_indices(input_pairs, sample_pairs, seed)
    selected_name_sha256 = write_sample(
        r1_path, r2_path, selected, input_pairs, r1_output, r2_output
    )
    manifest = {
        "algorithm": "sha256-counter-floyd-ordinal-v1",
        "cutadapt_report": Path(cutadapt_path).name,
        "input_fastq_1": Path(r1_path).name,
        "input_fastq_2": Path(r2_path).name,
        "input_pairs": input_pairs,
        "requested_pairs": sample_pairs,
        "seed": seed,
        "selected_name_sha256": selected_name_sha256,
        "selected_pairs": len(selected),
        "used_full_input": len(selected) == input_pairs,
    }
    Path(manifest_path).write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )


if __name__ == "__main__":
    main()
PYTHON

        screen_r1="~{fastqr1}"
        screen_r2="~{fastqr2}"
        if (( ~{sample_pairs} > 0 )); then
            screen_r1="~{SID}.sampled_R1.fastq.gz"
            screen_r2="~{SID}.sampled_R2.fastq.gz"
        fi

        mkdir -p genome/globin genome/rRNA genome/phix
        tar -xzf "~{globin_genome_dir_tar}" -C genome/globin --strip-components 1
        tar -xzf "~{rrna_genome_dir_tar}" -C genome/rRNA --strip-components 1
        tar -xzf "~{phix_genome_dir_tar}" -C genome/phix --strip-components 1

        run_screen() {
            kind="$1"
            genome_dir="$2"
            log="~{SID}_${kind}.log"
            report="~{SID}_${kind}_report.txt"

            bowtie2 \
                -p ~{ncpu} \
                -1 "${screen_r1}" \
                -2 "${screen_r2}" \
                -x "${genome_dir}/bowtie2_index" \
                --local \
                -S /dev/null \
                2> "${log}"

            awk -v id="~{SID}" -v kind="${kind}" '
                END {
                    if ($1 !~ /^[0-9]+([.][0-9]+)?%$/) {
                        print "invalid Bowtie2 terminal alignment percentage" > "/dev/stderr"
                        exit 2
                    }
                    print "Sample""\t""pct_"kind"\n"id"\t"$1
                }
            ' "${log}" > "${report}"
        }

        run_screen globin genome/globin
        run_screen rRNA genome/rRNA
        run_screen phix genome/phix
    >>>

    output {
        File globin_report = "${SID}_globin_report.txt"
        File rrna_report = "${SID}_rRNA_report.txt"
        File phix_report = "${SID}_phix_report.txt"
        Array[File] bowtie2_logs = [
            "${SID}_globin.log",
            "${SID}_rRNA.log",
            "${SID}_phix.log"
        ]
        File sampling_manifest = "${SID}_contamination_sampling.json"
    }

    runtime {
        cpu: ncpu
        memory: "${memory}GB"
        disks: "local-disk ${disk_space} HDD"
        docker: docker
        preemptible: preemptible
    }

    meta {
        description: "Run three contamination screens serially, optionally on a deterministic paired sample"
    }
}
