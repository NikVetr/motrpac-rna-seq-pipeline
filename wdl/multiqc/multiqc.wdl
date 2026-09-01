version 1.0

task multiQC {
    input {
        String SID
        Array[File] fastQCReports
        File trim_report

        Int memory
        Int disk_space
        Int ncpu
        Int preemptible
        String docker
    }

    command <<<
        set -euo pipefail
        echo "--- $(date "+[%b %d %H:%M:%S]") Beginning task, creating input directory ---"
        mkdir -p reports
        cd reports

        echo "--- $(date "+[%b %d %H:%M:%S]") Extracting fastQC reports from input tarball ---"
        for FILE in ~{sep=' ' fastQCReports}  ; do
            echo "Extracting $FILE"
            tar -xzf "$FILE"
        done

        cd ..

        echo "--- $(date "+[%b %d %H:%M:%S]") Running multiQC ---"
        multiqc \
          -d \
          -f \
          --cl-config "no_version_check: true" \
          -o multiQC_prealign_report \
          reports/* "~{trim_report}"

        echo "--- $(date "+[%b %d %H:%M:%S]") Creating output tarball ---"
        test -s multiQC_prealign_report/multiqc_report.html
        tar -czf "~{SID}.multiqc_prealign_report.tar.gz" ./multiQC_prealign_report

        echo "--- $(date "+[%b %d %H:%M:%S]") Finished creating output tarball, finished task ---"
    >>>

    output {
        File multiQC_report = "${SID}.multiqc_prealign_report.tar.gz"
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
            label: "Sample ID"
        }
        fastQCReports: {
            label: "FastQC reports"
        }
        trim_report: {
            label: "Trim report"
        }
    }
}
