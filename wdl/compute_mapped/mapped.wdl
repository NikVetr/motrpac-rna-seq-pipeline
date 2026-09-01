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
        Total=$(awk '{sum+=$3}END{print sum}' ~{SID}_aligned_chr_info.txt)
        grep "chrX" ~{SID}_aligned_chr_info.txt|awk -v tot="$Total"  -v name=~{SID} '{print "Sample""\t""pct_chrX""\n"name"\t"($3/tot)*100}' >chrX.txt
        grep "chrY" ~{SID}_aligned_chr_info.txt|awk -v tot="$Total" '{print "pct_chrY""\n"($3/tot)*100}' >chrY.txt
        grep "chrM" ~{SID}_aligned_chr_info.txt|awk -v tot="$Total" '{print "pct_chrM""\n"($3/tot)*100}' >chrM.txt
        grep "chr" ~{SID}_aligned_chr_info.txt|grep -v "chrX\|chrY\|chrM" |awk -v tot="$Total" '{sum+=$3}END{print "pct_chrAuto""\n"(sum/tot)*100}' >chrAuto.txt
        grep -v "chr\|^*" ~{SID}_aligned_chr_info.txt|awk -v tot="$Total" '{sum+=$3}END{print "pct_contig""\n"(sum/tot)*100}' >contig.txt

        echo "--- $(date "+[%b %d %H:%M:%S]") Consolidating intermediate files ---"
        paste chrX.txt chrY.txt chrM.txt chrAuto.txt contig.txt >"~{SID}_mapped_report.txt"

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
