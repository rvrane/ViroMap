#!/usr/bin/env python3

"""
remove_degenerate_bases.py

Remove degenerate bases from probe sequences by replacing them with G or C.
This script handles IUPAC degenerate nucleotide codes and replaces them
with the most GC-rich option to maintain probe stability.

Degenerate base replacements:
R (A/G) -> G
Y (C/T) -> C  
S (G/C) -> G (or C, both are GC)
W (A/T) -> C (choosing GC-rich option)
K (G/T) -> G
M (A/C) -> C
B (C/G/T) -> G
D (A/G/T) -> G
H (A/C/T) -> C
V (A/C/G) -> G
N (A/C/G/T) -> G
"""

import os
import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import glob

def replace_degenerate_bases(sequence):
    """
    Replace degenerate bases with G or C to maintain high GC content.
    
    Args:
        sequence (str): DNA sequence with potential degenerate bases
        
    Returns:
        str: Sequence with degenerate bases replaced
    """
    
    # Define replacement mapping - prioritizing G and C
    degenerate_map = {
        'R': 'G',  # R = A/G -> choose G
        'Y': 'C',  # Y = C/T -> choose C
        'S': 'G',  # S = G/C -> choose G (could also be C)
        'W': 'C',  # W = A/T -> choose C (GC-rich option)
        'K': 'G',  # K = G/T -> choose G
        'M': 'C',  # M = A/C -> choose C
        'B': 'G',  # B = C/G/T -> choose G
        'D': 'G',  # D = A/G/T -> choose G
        'H': 'C',  # H = A/C/T -> choose C
        'V': 'G',  # V = A/C/G -> choose G
        'N': 'G',  # N = A/C/G/T -> choose G
    }
    
    # Convert to uppercase and replace degenerate bases
    clean_sequence = sequence.upper()
    original_length = len(clean_sequence)
    
    # Count degenerate bases before replacement
    degenerate_count = {}
    for base in degenerate_map.keys():
        count = clean_sequence.count(base)
        if count > 0:
            degenerate_count[base] = count
    
    # Replace degenerate bases
    for degenerate, replacement in degenerate_map.items():
        clean_sequence = clean_sequence.replace(degenerate, replacement)
    
    return clean_sequence, degenerate_count, original_length

def calculate_gc_content(sequence):
    """Calculate GC content percentage."""
    if len(sequence) == 0:
        return 0
    gc_count = sequence.count('G') + sequence.count('C')
    return (gc_count / len(sequence)) * 100

def process_fasta_file(input_fasta, output_fasta):
    """
    Process a single FASTA file to remove degenerate bases.
    
    Args:
        input_fasta (str): Input FASTA file path
        output_fasta (str): Output FASTA file path
        
    Returns:
        dict: Statistics about the processing
    """
    
    if not os.path.exists(input_fasta):
        print(f"Error: Input file not found: {input_fasta}")
        return None
    
    processed_records = []
    stats = {
        'total_probes': 0,
        'probes_with_degenerates': 0,
        'total_degenerate_bases': 0,
        'degenerate_types': {},
        'gc_content_before': [],
        'gc_content_after': []
    }
    
    print(f"Processing: {os.path.basename(input_fasta)}")
    
    try:
        for record in SeqIO.parse(input_fasta, "fasta"):
            stats['total_probes'] += 1
            original_seq = str(record.seq).upper()
            
            # Calculate original GC content
            original_gc = calculate_gc_content(original_seq)
            stats['gc_content_before'].append(original_gc)
            
            # Replace degenerate bases
            clean_seq, degenerate_count, original_length = replace_degenerate_bases(original_seq)
            
            # Calculate new GC content
            new_gc = calculate_gc_content(clean_seq)
            stats['gc_content_after'].append(new_gc)
            
            # Update statistics
            if degenerate_count:
                stats['probes_with_degenerates'] += 1
                for base, count in degenerate_count.items():
                    stats['total_degenerate_bases'] += count
                    if base in stats['degenerate_types']:
                        stats['degenerate_types'][base] += count
                    else:
                        stats['degenerate_types'][base] = count
            
            # Create new record
            new_record = SeqRecord(
                Seq(clean_seq),
                id=record.id.replace('_with_adapters', '_clean'),
                description=f"{record.description} | Degenerates_removed | GC: {new_gc:.1f}%"
            )
            
            processed_records.append(new_record)
            
            # Show details for first few probes with degenerates
            if degenerate_count and stats['probes_with_degenerates'] <= 3:
                print(f"  {record.id}:")
                print(f"    Degenerates found: {degenerate_count}")
                print(f"    GC content: {original_gc:.1f}% -> {new_gc:.1f}%")
        
        # Write processed sequences
        with open(output_fasta, 'w') as output_handle:
            SeqIO.write(processed_records, output_handle, "fasta")
        
        return stats
        
    except Exception as e:
        print(f"Error processing {input_fasta}: {str(e)}")
        return None

def process_directory(input_dir, output_dir=None):
    """
    Process all FASTA files in a directory.
    
    Args:
        input_dir (str): Input directory containing FASTA files
        output_dir (str): Output directory (if None, creates 'clean_probes' subdirectory)
    """
    
    if not os.path.exists(input_dir):
        print(f"Error: Directory not found: {input_dir}")
        return
    
    # Set output directory
    if output_dir is None:
        output_dir = os.path.join(input_dir, "clean_probes")
    
    os.makedirs(output_dir, exist_ok=True)
    
    # Find all FASTA files
    fasta_patterns = ['*.fasta', '*.fa', '*.fas']
    fasta_files = []
    for pattern in fasta_patterns:
        fasta_files.extend(glob.glob(os.path.join(input_dir, pattern)))
    
    if not fasta_files:
        print(f"No FASTA files found in {input_dir}")
        return
    
    print(f"Found {len(fasta_files)} FASTA files to process")
    print(f"Output directory: {output_dir}")
    print("="*70)
    
    # Process all files
    total_stats = {
        'total_files': len(fasta_files),
        'successful_files': 0,
        'total_probes': 0,
        'probes_with_degenerates': 0,
        'total_degenerate_bases': 0,
        'degenerate_types': {},
        'avg_gc_before': 0,
        'avg_gc_after': 0
    }
    
    all_gc_before = []
    all_gc_after = []
    
    for fasta_file in sorted(fasta_files):
        input_path = fasta_file
        base_name = os.path.splitext(os.path.basename(fasta_file))[0]
        output_filename = f"{base_name}_clean.fasta"
        output_path = os.path.join(output_dir, output_filename)
        
        # Process file
        stats = process_fasta_file(input_path, output_path)
        
        if stats:
            total_stats['successful_files'] += 1
            total_stats['total_probes'] += stats['total_probes']
            total_stats['probes_with_degenerates'] += stats['probes_with_degenerates']
            total_stats['total_degenerate_bases'] += stats['total_degenerate_bases']
            
            # Merge degenerate types
            for base, count in stats['degenerate_types'].items():
                if base in total_stats['degenerate_types']:
                    total_stats['degenerate_types'][base] += count
                else:
                    total_stats['degenerate_types'][base] = count
            
            all_gc_before.extend(stats['gc_content_before'])
            all_gc_after.extend(stats['gc_content_after'])
            
            print(f"  ✅ Processed: {stats['total_probes']} probes, "
                  f"{stats['probes_with_degenerates']} had degenerates")
        else:
            print(f"  ❌ Failed to process: {os.path.basename(fasta_file)}")
        
        print("-"*50)
    
    # Calculate averages
    if all_gc_before:
        total_stats['avg_gc_before'] = sum(all_gc_before) / len(all_gc_before)
    if all_gc_after:
        total_stats['avg_gc_after'] = sum(all_gc_after) / len(all_gc_after)
    
    # Print summary
    print("\n" + "="*70)
    print("DEGENERATE BASE REMOVAL SUMMARY")
    print("="*70)
    print(f"Files processed successfully: {total_stats['successful_files']}/{total_stats['total_files']}")
    print(f"Total probes processed: {total_stats['total_probes']}")
    print(f"Probes with degenerate bases: {total_stats['probes_with_degenerates']}")
    print(f"Total degenerate bases removed: {total_stats['total_degenerate_bases']}")
    print(f"Average GC content: {total_stats['avg_gc_before']:.1f}% -> {total_stats['avg_gc_after']:.1f}%")
    
    if total_stats['degenerate_types']:
        print(f"\nDegenerate bases found and replaced:")
        for base, count in sorted(total_stats['degenerate_types'].items()):
            replacement = {'R': 'G', 'Y': 'C', 'S': 'G', 'W': 'C', 'K': 'G', 
                          'M': 'C', 'B': 'G', 'D': 'G', 'H': 'C', 'V': 'G', 'N': 'G'}
            print(f"  {base} -> {replacement.get(base, '?')}: {count:,} bases")
    
    print(f"\nCleaned probe files saved in: {output_dir}")
    print("="*70)

def main():
    if len(sys.argv) < 2:
        print("Usage:")
        print("  Single file: python3 remove_degenerate_bases.py <input.fasta> [output.fasta]")
        print("  Directory:   python3 remove_degenerate_bases.py <probe_directory> [output_directory]")
        print()
        print("Examples:")
        print("  python3 remove_degenerate_bases.py probes.fasta probes_clean.fasta")
        print("  python3 remove_degenerate_bases.py /path/to/probe_results/")
        print()
        print("Degenerate bases will be replaced as follows:")
        print("  R (A/G) -> G,  Y (C/T) -> C,  S (G/C) -> G")
        print("  W (A/T) -> C,  K (G/T) -> G,  M (A/C) -> C")
        print("  B (C/G/T) -> G,  D (A/G/T) -> G,  H (A/C/T) -> C")
        print("  V (A/C/G) -> G,  N (A/C/G/T) -> G")
        sys.exit(1)
    
    input_path = sys.argv[1]
    
    # Check if input is a file or directory
    if os.path.isfile(input_path):
        # Single file mode
        if len(sys.argv) >= 3:
            output_path = sys.argv[2]
        else:
            # Generate output filename
            base_name = os.path.splitext(input_path)[0]
            output_path = f"{base_name}_clean.fasta"
        
        print("Degenerate Base Removal Tool")
        print("="*40)
        
        stats = process_fasta_file(input_path, output_path)
        if stats:
            print(f"\n✅ Success! Clean probes saved to: {output_path}")
            if stats['probes_with_degenerates'] > 0:
                print(f"Removed {stats['total_degenerate_bases']} degenerate bases from {stats['probes_with_degenerates']} probes")
            else:
                print("No degenerate bases found - all probes were already clean")
        
    elif os.path.isdir(input_path):
        # Directory mode
        output_dir = sys.argv[2] if len(sys.argv) >= 3 else None
        
        print("Degenerate Base Removal Tool - Directory Mode")
        print("="*50)
        
        process_directory(input_path, output_dir)
        
    else:
        print(f"Error: Path not found: {input_path}")
        sys.exit(1)

if __name__ == "__main__":
    main()