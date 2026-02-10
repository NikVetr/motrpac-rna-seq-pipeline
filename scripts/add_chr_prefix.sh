#!/usr/bin/env bash
# Add UCSC-style 'chr' prefix to Ensembl chromosome names in FASTA and GTF files.
#
# The pipeline expects chr-prefixed chromosome names (chr1, chrX, chrM, etc.)
# but Ensembl uses bare names (1, X, MT). This script converts Ensembl naming
# to UCSC/chr convention to match the format used for rn6 and rn7.
#
# Mapping:
#   1, 2, ..., 20, X, Y  →  chr1, chr2, ..., chr20, chrX, chrY
#   MT                    →  chrM
#   Scaffolds/contigs     →  unchanged
#
# Usage: bash scripts/add_chr_prefix.sh <fasta_file> <gtf_file>

set -euo pipefail

if [ $# -ne 2 ]; then
    echo "Usage: $0 <fasta_file> <gtf_file>"
    echo "Example: $0 data/rat-ensembl-release-115/Rattus_norvegicus.GRCr8.dna.toplevel.fa data/rat-ensembl-release-115/Rattus_norvegicus.GRCr8.115.gtf"
    exit 1
fi

FASTA="$1"
GTF="$2"

# Verify input files exist
for f in "$FASTA" "$GTF"; do
    if [ ! -f "$f" ]; then
        echo "Error: File not found: $f"
        exit 1
    fi
done

echo "Adding chr prefix to FASTA: $FASTA"
sed -E 's/^>([0-9]+|X|Y) />chr\1 /; s/^>MT />chrM /' "$FASTA" > "${FASTA}.tmp" && mv "${FASTA}.tmp" "$FASTA"

echo "Adding chr prefix to GTF: $GTF"
sed -E '/^#/! { s/^([0-9]+|X|Y)\t/chr\1\t/; s/^MT\t/chrM\t/; }' "$GTF" > "${GTF}.tmp" && mv "${GTF}.tmp" "$GTF"

echo ""
echo "Done. Verify chromosome names:"
echo ""
echo "FASTA headers (first 5):"
grep '^>' "$FASTA" | head -5
echo ""
echo "GTF chromosomes (unique):"
grep -v '^#' "$GTF" | cut -f1 | sort -u | head -10
