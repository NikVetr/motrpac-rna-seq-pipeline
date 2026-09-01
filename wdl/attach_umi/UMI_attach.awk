#!/usr/bin/awk -f

# Append an eight-base index-read UMI while requiring complete, synchronized
# biological and index FASTQ records. Valid output matches the legacy format.

BEGIN {
    FS = " "
    if (Ifq == "") {
        Ifq = "index.fastq.gz"
    }
    index_command = "gzip -cd -- " Ifq
    failed = 0
    records = 0
}

function fail(message) {
    print "UMI_attach: " message > "/dev/stderr"
    failed = 1
    exit 1
}

function next_index_line(context, status) {
    status = (index_command | getline index_line)
    if (status == 0) {
        fail("index FASTQ ended while reading " context " for biological record " (records + 1))
    }
    if (status < 0) {
        fail("could not read " context " from index FASTQ")
    }
    return index_line
}

NR % 4 == 1 {
    if (substr($0, 1, 1) != "@") {
        fail("biological FASTQ record " (records + 1) " has an invalid header")
    }

    read_header = $0
    read_name = $1
    index_header = next_index_line("header")
    if (substr(index_header, 1, 1) != "@") {
        fail("index FASTQ record " (records + 1) " has an invalid header")
    }
    split(index_header, index_header_fields, " ")
    if (index_header_fields[1] != read_name) {
        fail("record-name mismatch at biological FASTQ line " NR ": " read_name " != " index_header_fields[1])
    }

    umi = next_index_line("sequence")
    index_plus = next_index_line("separator")
    index_quality = next_index_line("quality")
    if (substr(index_plus, 1, 1) != "+") {
        fail("index FASTQ record " (records + 1) " has an invalid separator")
    }
    if (length(umi) == 0 || length(umi) != length(index_quality)) {
        fail("index FASTQ record " (records + 1) " has inconsistent sequence and quality lengths")
    }
    if (length(umi) != 8) {
        fail("index FASTQ record " (records + 1) " has UMI length " length(umi) "; expected 8")
    }

    print read_name ":" umi substr(read_header, length(read_name) + 1)
}

NR % 4 == 2 {
    read_sequence = $0
    print
}

NR % 4 == 3 {
    if (substr($0, 1, 1) != "+") {
        fail("biological FASTQ record " (records + 1) " has an invalid separator")
    }
    print "+"
}

NR % 4 == 0 {
    if (length(read_sequence) != length($0)) {
        fail("biological FASTQ record " (records + 1) " has inconsistent sequence and quality lengths")
    }
    print
    records++
}

END {
    if (failed) {
        close(index_command)
        exit 1
    }
    if (NR == 0 || NR % 4 != 0) {
        close(index_command)
        print "UMI_attach: biological FASTQ is empty or truncated" > "/dev/stderr"
        exit 1
    }

    extra_status = (index_command | getline extra_index_line)
    if (extra_status > 0) {
        close(index_command)
        print "UMI_attach: index FASTQ contains more records than the biological FASTQ" > "/dev/stderr"
        exit 1
    }
    if (extra_status < 0) {
        close(index_command)
        print "UMI_attach: failed while checking the end of the index FASTQ" > "/dev/stderr"
        exit 1
    }
    if (close(index_command) != 0) {
        print "UMI_attach: index FASTQ decompression or checksum validation failed" > "/dev/stderr"
        exit 1
    }
}
