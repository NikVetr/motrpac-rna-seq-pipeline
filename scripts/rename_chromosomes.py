#!/usr/bin/env python3
"""
Rename chromosome names from NCBI RefSeq accessions to UCSC-style names.

This script converts NCBI RefSeq chromosome naming (e.g., NC_086019.1) to 
UCSC-style naming (e.g., chr1) using the sequence_report.jsonl file provided
by NCBI Datasets downloads.

Usage:
    python rename_chromosomes.py <sequence_report.jsonl> <input_dir> <output_dir>

Example:
    python rename_chromosomes.py \
        data/GCF_036323735.1/sequence_report.jsonl \
        data/GCF_036323735.1 \
        data/GCF_036323735.1_ucsc

The script will:
1. Generate a chromosome name mapping from sequence_report.jsonl
2. Rename chromosomes in the genomic FASTA file
3. Rename chromosomes in the GTF annotation file
4. Perform QC checks to verify the renaming was successful
5. Generate a summary plot of genes per chromosome
"""

import argparse
import json
import gzip
import sys
import os
import re
from pathlib import Path
from collections import defaultdict


def load_chromosome_mapping(sequence_report_path: str) -> dict:
    """
    Load chromosome name mapping from NCBI sequence_report.jsonl.
    
    Returns a dict mapping RefSeq accession -> UCSC-style name
    """
    mapping = {}
    
    with open(sequence_report_path, 'r') as f:
        for line in f:
            record = json.loads(line)
            refseq_acc = record.get('refseqAccession')
            ucsc_name = record.get('sequenceName')
            
            if refseq_acc and ucsc_name:
                mapping[refseq_acc] = ucsc_name
    
    return mapping


def open_file(filepath: str, mode: str = 'r'):
    """Open a file, handling gzip compression automatically."""
    if filepath.endswith('.gz'):
        return gzip.open(filepath, mode + 't')
    return open(filepath, mode)


def rename_fasta(input_fasta: str, output_fasta: str, mapping: dict) -> dict:
    """
    Rename chromosome names in a FASTA file.
    
    Returns stats about the renaming operation.
    """
    stats = {
        'total_sequences': 0,
        'renamed_sequences': 0,
        'unmapped_sequences': [],
        'input_size': os.path.getsize(input_fasta),
    }
    
    with open_file(input_fasta, 'r') as fin, open(output_fasta, 'w') as fout:
        for line in fin:
            if line.startswith('>'):
                stats['total_sequences'] += 1
                # Extract the accession from the header (first word after >)
                header_parts = line[1:].split()
                accession = header_parts[0]
                
                if accession in mapping:
                    new_name = mapping[accession]
                    # Preserve the rest of the header as a comment
                    rest_of_header = ' '.join(header_parts[1:]) if len(header_parts) > 1 else ''
                    if rest_of_header:
                        new_header = f">{new_name} {rest_of_header}\n"
                    else:
                        new_header = f">{new_name}\n"
                    fout.write(new_header)
                    stats['renamed_sequences'] += 1
                else:
                    # Keep original if no mapping exists
                    fout.write(line)
                    stats['unmapped_sequences'].append(accession)
            else:
                fout.write(line)
    
    stats['output_size'] = os.path.getsize(output_fasta)
    return stats


def rename_gtf(input_gtf: str, output_gtf: str, mapping: dict) -> dict:
    """
    Rename chromosome names in a GTF file (first column).
    
    Returns stats about the renaming operation.
    """
    stats = {
        'total_lines': 0,
        'renamed_lines': 0,
        'unmapped_chromosomes': set(),
        'genes_per_chromosome': defaultdict(set),
    }
    
    with open_file(input_gtf, 'r') as fin, open(output_gtf, 'w') as fout:
        for line in fin:
            stats['total_lines'] += 1
            
            # Skip comment lines
            if line.startswith('#'):
                fout.write(line)
                continue
            
            parts = line.split('\t')
            if len(parts) >= 9:
                old_chrom = parts[0]
                
                if old_chrom in mapping:
                    new_chrom = mapping[old_chrom]
                    parts[0] = new_chrom
                    fout.write('\t'.join(parts))
                    stats['renamed_lines'] += 1
                    
                    # Track genes per chromosome for QC plot
                    feature_type = parts[2]
                    if feature_type in ('gene', 'transcript'):
                        # Extract gene_id from attributes
                        attrs = parts[8]
                        gene_match = re.search(r'gene_id\s+"([^"]+)"', attrs)
                        if gene_match:
                            stats['genes_per_chromosome'][new_chrom].add(gene_match.group(1))
                else:
                    fout.write(line)
                    stats['unmapped_chromosomes'].add(old_chrom)
            else:
                fout.write(line)
    
    # Convert sets to counts
    stats['genes_per_chromosome'] = {k: len(v) for k, v in stats['genes_per_chromosome'].items()}
    stats['unmapped_chromosomes'] = list(stats['unmapped_chromosomes'])
    
    return stats


def generate_qc_report(fasta_stats: dict, gtf_stats: dict, output_dir: str) -> str:
    """Generate a QC report as text."""
    report = []
    report.append("=" * 60)
    report.append("CHROMOSOME RENAMING QC REPORT")
    report.append("=" * 60)
    report.append("")
    
    # FASTA QC
    report.append("FASTA File:")
    report.append(f"  Total sequences: {fasta_stats['total_sequences']}")
    report.append(f"  Renamed sequences: {fasta_stats['renamed_sequences']}")
    report.append(f"  Input file size: {fasta_stats['input_size']:,} bytes")
    report.append(f"  Output file size: {fasta_stats['output_size']:,} bytes")
    
    # Size check - should be very similar (small differences due to header changes)
    size_diff_pct = abs(fasta_stats['output_size'] - fasta_stats['input_size']) / fasta_stats['input_size'] * 100
    if size_diff_pct < 1:
        report.append(f"  ✓ File size difference: {size_diff_pct:.3f}% (OK)")
    else:
        report.append(f"  ✗ WARNING: File size difference: {size_diff_pct:.2f}%")
    
    if fasta_stats['unmapped_sequences']:
        report.append(f"  ⚠ Unmapped sequences ({len(fasta_stats['unmapped_sequences'])}): {', '.join(fasta_stats['unmapped_sequences'][:5])}...")
    else:
        report.append("  ✓ All sequences were mapped")
    
    report.append("")
    
    # GTF QC
    report.append("GTF File:")
    report.append(f"  Total lines: {gtf_stats['total_lines']:,}")
    report.append(f"  Renamed lines: {gtf_stats['renamed_lines']:,}")
    
    if gtf_stats['unmapped_chromosomes']:
        report.append(f"  ⚠ Unmapped chromosomes: {', '.join(gtf_stats['unmapped_chromosomes'][:5])}")
    else:
        report.append("  ✓ All chromosomes were mapped")
    
    report.append("")
    
    # Gene counts per chromosome
    report.append("Genes per Chromosome:")
    genes_per_chr = gtf_stats['genes_per_chromosome']
    
    # Sort chromosomes naturally (chr1, chr2, ..., chr10, ..., chrX, chrY)
    def chr_sort_key(x):
        if x.startswith('chr'):
            rest = x[3:]
            if rest.isdigit():
                return (0, int(rest))
            elif rest == 'X':
                return (1, 0)
            elif rest == 'Y':
                return (1, 1)
            elif rest == 'M':
                return (1, 2)
            else:
                return (2, rest)
        return (3, x)
    
    sorted_chroms = sorted(genes_per_chr.keys(), key=chr_sort_key)
    
    for chrom in sorted_chroms:
        count = genes_per_chr[chrom]
        bar = '█' * min(count // 500, 50)  # Scale bar
        report.append(f"  {chrom:25s} {count:6d} {bar}")
    
    report.append("")
    report.append(f"Total genes: {sum(genes_per_chr.values()):,}")
    report.append("=" * 60)
    
    return '\n'.join(report)


def generate_qc_plot(gtf_stats: dict, output_path: str):
    """Generate a bar plot of genes per chromosome."""
    try:
        import matplotlib
        matplotlib.use('Agg')  # Non-interactive backend
        import matplotlib.pyplot as plt
    except ImportError:
        print("Warning: matplotlib not available, skipping plot generation")
        return
    
    genes_per_chr = gtf_stats['genes_per_chromosome']
    
    # Sort chromosomes naturally
    def chr_sort_key(x):
        if x.startswith('chr'):
            rest = x[3:]
            if rest.isdigit():
                return (0, int(rest))
            elif rest == 'X':
                return (1, 0)
            elif rest == 'Y':
                return (1, 1)
            elif rest == 'M':
                return (1, 2)
            else:
                return (2, rest)
        return (3, x)
    
    # Filter to main chromosomes for cleaner plot
    main_chroms = [c for c in genes_per_chr.keys() 
                   if re.match(r'^chr(\d+|X|Y|M)$', c)]
    main_chroms = sorted(main_chroms, key=chr_sort_key)
    
    counts = [genes_per_chr[c] for c in main_chroms]
    
    # Create figure
    fig, ax = plt.subplots(figsize=(14, 6))
    
    bars = ax.bar(range(len(main_chroms)), counts, color='steelblue', edgecolor='navy')
    
    ax.set_xticks(range(len(main_chroms)))
    ax.set_xticklabels([c.replace('chr', '') for c in main_chroms], fontsize=10)
    ax.set_xlabel('Chromosome', fontsize=12)
    ax.set_ylabel('Number of Genes', fontsize=12)
    ax.set_title('Gene Count per Chromosome (after UCSC-style renaming)', fontsize=14)
    
    # Add value labels on bars
    for bar, count in zip(bars, counts):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 50, 
                f'{count:,}', ha='center', va='bottom', fontsize=8, rotation=90)
    
    ax.set_ylim(0, max(counts) * 1.15)
    ax.grid(axis='y', alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"QC plot saved to: {output_path}")


def save_mapping_file(mapping: dict, output_path: str):
    """Save the chromosome mapping to a TSV file for reference."""
    with open(output_path, 'w') as f:
        f.write("# RefSeq_Accession\tUCSC_Name\n")
        for refseq, ucsc in sorted(mapping.items()):
            f.write(f"{refseq}\t{ucsc}\n")


def main():
    parser = argparse.ArgumentParser(
        description='Rename chromosome names from NCBI RefSeq to UCSC-style naming',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__
    )
    parser.add_argument('sequence_report', 
                        help='Path to sequence_report.jsonl from NCBI Datasets')
    parser.add_argument('input_dir',
                        help='Input directory containing FASTA and GTF files')
    parser.add_argument('output_dir',
                        help='Output directory for renamed files')
    parser.add_argument('--fasta-pattern', default='*_genomic.fna',
                        help='Glob pattern for genomic FASTA file (default: *_genomic.fna)')
    parser.add_argument('--gtf-pattern', default='*.gtf',
                        help='Glob pattern for GTF file (default: *.gtf)')
    parser.add_argument('--skip-plot', action='store_true',
                        help='Skip generating the QC plot')
    parser.add_argument('--add-suffix', action='store_true',
                        help='Add _ucsc suffix to output filenames (default: keep original names)')
    
    args = parser.parse_args()
    
    # Resolve paths
    input_dir = Path(args.input_dir)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Find input files
    fasta_files = list(input_dir.glob(args.fasta_pattern))
    gtf_files = list(input_dir.glob(args.gtf_pattern))
    
    if not fasta_files:
        print(f"Error: No FASTA file matching '{args.fasta_pattern}' found in {input_dir}")
        sys.exit(1)
    if not gtf_files:
        print(f"Error: No GTF file matching '{args.gtf_pattern}' found in {input_dir}")
        sys.exit(1)
    
    input_fasta = fasta_files[0]
    input_gtf = gtf_files[0]
    
    print(f"Input FASTA: {input_fasta}")
    print(f"Input GTF: {input_gtf}")
    print(f"Output directory: {output_dir}")
    print()
    
    # Load chromosome mapping
    print("Loading chromosome mapping from sequence_report.jsonl...")
    mapping = load_chromosome_mapping(args.sequence_report)
    print(f"  Loaded {len(mapping)} chromosome mappings")
    
    # Save mapping file for reference
    mapping_file = output_dir / 'chromosome_mapping.tsv'
    save_mapping_file(mapping, mapping_file)
    print(f"  Mapping saved to: {mapping_file}")
    print()
    
    # Determine output filenames
    if args.add_suffix:
        output_fasta = output_dir / input_fasta.name.replace('.fna', '_ucsc.fna')
        output_gtf = output_dir / input_gtf.name.replace('.gtf', '_ucsc.gtf')
    else:
        output_fasta = output_dir / input_fasta.name
        output_gtf = output_dir / input_gtf.name
    
    # Rename FASTA
    print(f"Renaming FASTA chromosomes -> {output_fasta.name}...")
    fasta_stats = rename_fasta(str(input_fasta), str(output_fasta), mapping)
    print(f"  Processed {fasta_stats['total_sequences']} sequences")
    print()
    
    # Rename GTF
    print(f"Renaming GTF chromosomes -> {output_gtf.name}...")
    gtf_stats = rename_gtf(str(input_gtf), str(output_gtf), mapping)
    print(f"  Processed {gtf_stats['total_lines']:,} lines")
    print()
    
    # Generate QC report
    qc_report = generate_qc_report(fasta_stats, gtf_stats, str(output_dir))
    print(qc_report)
    
    # Save QC report
    qc_report_path = output_dir / 'qc_report.txt'
    with open(qc_report_path, 'w') as f:
        f.write(qc_report)
    print(f"\nQC report saved to: {qc_report_path}")
    
    # Generate QC plot
    if not args.skip_plot:
        plot_path = output_dir / 'genes_per_chromosome.png'
        generate_qc_plot(gtf_stats, str(plot_path))
    
    # Final validation
    print("\n" + "=" * 60)
    print("VALIDATION SUMMARY")
    print("=" * 60)
    
    all_ok = True
    
    # Check sequence count matches
    if fasta_stats['total_sequences'] == fasta_stats['renamed_sequences'] + len(fasta_stats['unmapped_sequences']):
        print("✓ FASTA sequence count is consistent")
    else:
        print("✗ FASTA sequence count mismatch!")
        all_ok = False
    
    # Check file sizes are similar
    size_diff_pct = abs(fasta_stats['output_size'] - fasta_stats['input_size']) / fasta_stats['input_size'] * 100
    if size_diff_pct < 1:
        print(f"✓ FASTA file sizes match within 1% ({size_diff_pct:.3f}%)")
    else:
        print(f"✗ FASTA file size differs by {size_diff_pct:.2f}%")
        all_ok = False
    
    # Check main chromosomes are present
    expected_main_chroms = {'chr1', 'chr2', 'chrX', 'chrY'}
    found_chroms = set(gtf_stats['genes_per_chromosome'].keys())
    if expected_main_chroms.issubset(found_chroms):
        print("✓ Main chromosomes (chr1, chr2, chrX, chrY) found in output")
    else:
        missing = expected_main_chroms - found_chroms
        print(f"✗ Missing expected chromosomes: {missing}")
        all_ok = False
    
    if all_ok:
        print("\n✓ All validation checks passed!")
        print(f"\nOutput files ready in: {output_dir}/")
        print(f"  - {output_fasta.name}")
        print(f"  - {output_gtf.name}")
        print(f"  - chromosome_mapping.tsv")
        print(f"  - qc_report.txt")
        if not args.skip_plot:
            print(f"  - genes_per_chromosome.png")
    else:
        print("\n⚠ Some validation checks failed - please review the output")
        sys.exit(1)


if __name__ == '__main__':
    main()
