version 1.0

task samtools_mapped {
    input {
        String SID
        File input_bam

        Int memory
        Int disk_space
        Int ncpu
        Int preemptible
        String docker
    }

    command <<<
        set -euo pipefail
        echo "--- $(date "+[%b %d %H:%M:%S]") Beginning task ---"

        echo "--- $(date "+[%b %d %H:%M:%S]") Counting primary alignments by reference ---"
        samtools view -H ~{input_bam} > ~{SID}_header.sam
        samtools view -F 0x900 ~{input_bam} | awk '
            BEGIN {
                FS = OFS = "\t"
                while ((getline line < "~{SID}_header.sam") > 0) {
                    field_count = split(line, fields, "\t")
                    if (fields[1] == "@SQ") {
                        name = ""
                        reference_size = ""
                        for (i = 2; i <= field_count; i++) {
                            if (fields[i] ~ /^SN:/) {
                                name = substr(fields[i], 4)
                            } else if (fields[i] ~ /^LN:/) {
                                reference_size = substr(fields[i], 4)
                            }
                        }
                        if (name == "" || reference_size == "") {
                            print "invalid @SQ header line: " line > "/dev/stderr"
                            exit 2
                        }
                        key = "R:" name
                        if (key in known_reference) {
                            print "duplicate @SQ reference: " name > "/dev/stderr"
                            exit 2
                        }
                        reference_order[++reference_count] = name
                        reference_length[key] = reference_size
                        known_reference[key] = 1
                    }
                }
                close("~{SID}_header.sam")
            }
            {
                flag = $2 + 0
                reference = $3
                key = (reference == "*" ? "*" : "R:" reference)
                if (reference != "*" && !(key in known_reference)) {
                    print "alignment has unknown reference: " reference > "/dev/stderr"
                    exit 2
                }
                if (int(flag / 4) % 2 == 1) {
                    unmapped[key]++
                } else {
                    mapped[key]++
                }
            }
            END {
                for (i = 1; i <= reference_count; i++) {
                    reference = reference_order[i]
                    key = "R:" reference
                    printf "%s\t%s\t%.0f\t%.0f\n", reference, reference_length[key], mapped[key] + 0, unmapped[key] + 0
                }
                printf "*\t0\t0\t%.0f\n", unmapped["*"] + 0
            }
        ' > ~{SID}_aligned_chr_info.txt

        echo "--- $(date "+[%b %d %H:%M:%S]") Extracting reports, info ---"
        awk -v name=~{SID} '
            BEGIN { FS = OFS = "\t" }
            $1 != "*" {
                count = $3 + 0
                total += count
                if ($1 == "chrX") {
                    chr_x += count
                } else if ($1 == "chrY") {
                    chr_y += count
                } else if ($1 == "chrM") {
                    chr_m += count
                } else if ($1 ~ /^chr[0-9]+$/) {
                    chr_auto += count
                } else {
                    contig += count
                }
            }
            END {
                if (total == 0) {
                    print "no mapped primary alignments" > "/dev/stderr"
                    exit 2
                }
                print "Sample", "pct_chrX", "pct_chrY", "pct_chrM", "pct_chrAuto", "pct_contig"
                print name, chr_x / total * 100, chr_y / total * 100, chr_m / total * 100, chr_auto / total * 100, contig / total * 100
            }
        ' ~{SID}_aligned_chr_info.txt > "~{SID}_mapped_report.txt"

        echo "--- $(date "+[%b %d %H:%M:%S]") Task complete ---"
    >>>

    output {
        File aligned_chrinfo = "${SID}_aligned_chr_info.txt"
        File report = "${SID}_mapped_report.txt"
    }

    runtime {
        cpu: ncpu
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
            label: "Aligned BAM File"
        }
    }

    meta {
        author: "Archana Raja"
    }
}
