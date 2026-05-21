import os
from Bio import SeqIO
import pandas as pd

def calc_gc_content(seq):
    gc = seq.count('G') + seq.count('C')
    return 100.0 * gc / len(seq) if len(seq) > 0 else 0

def melting_temp(seq):
    at = seq.count('A') + seq.count('T')
    gc = seq.count('G') + seq.count('C')
    return 2 * at + 4 * gc

def has_homopolymer(seq, max_len=7):
    count, prev = 1, ''
    for base in seq:
        if base == prev:
            count += 1
            if count > max_len:
                return True
        else:
            count = 1
            prev = base
    return False

def design_probes_consecutive(input_fasta, output_dir, probe_length=120):  # Changed from 80 to 120
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    probe_records = []
    summary_rows = []

    # Set permissive filters
    min_gc = 20
    max_gc = 80
    max_homopolymer = 7

    for record in SeqIO.parse(input_fasta, "fasta"):
        seq = str(record.seq).upper()
        seq_len = len(seq)
        probes = []
        reasons_skipped = {"ambig": 0, "gc": 0, "homopolymer": 0}
        total_probes_attempted = 0

        if seq_len < probe_length:
            print(f"Sequence {record.id} too short: {seq_len}bp, skipping.")
            summary_rows.append({
                'sequence_id': record.id,
                'num_probes': 0,
                'mean_gc_content': None,
                'mean_tm': None
            })
            continue

        position = 0
        while position < seq_len:
            # Determine probe end position
            end_pos = min(position + probe_length, seq_len)
            current_length = end_pos - position

            # If we have less than 120bp left, borrow from the left
            if current_length < probe_length and position > 0:
                needed = probe_length - current_length
                start_pos = max(0, position - needed)
                probe_seq = seq[start_pos:end_pos]
                if len(probe_seq) > probe_length:
                    probe_seq = probe_seq[-probe_length:]
            elif current_length >= probe_length:
                probe_seq = seq[position:position + probe_length]
            else:
                probe_seq = seq[position:end_pos]
                if len(probe_seq) < probe_length:
                    print(f"Remaining sequence too short for probe: {len(probe_seq)}bp")
                    break

            total_probes_attempted += 1

            # Apply filters
            if 'N' in probe_seq or 'X' in probe_seq:
                reasons_skipped["ambig"] += 1
                position += probe_length
                continue

            gc = calc_gc_content(probe_seq)
            if gc < min_gc or gc > max_gc:
                reasons_skipped["gc"] += 1
                position += probe_length
                continue

            if has_homopolymer(probe_seq, max_len=max_homopolymer):
                reasons_skipped["homopolymer"] += 1
                position += probe_length
                continue

            probes.append(probe_seq)
            probe_records.append((record.id, probe_seq, position))

            position += probe_length

        if probes:
            gc_values = [calc_gc_content(p) for p in probes]
            tm_values = [melting_temp(p) for p in probes]
            summary_rows.append({
                'sequence_id': record.id,
                'num_probes': len(probes),
                'mean_gc_content': sum(gc_values) / len(gc_values),
                'mean_tm': sum(tm_values) / len(tm_values)
            })
        else:
            summary_rows.append({
                'sequence_id': record.id,
                'num_probes': 0,
                'mean_gc_content': None,
                'mean_tm': None
            })

        print(f"{record.id}: attempted={total_probes_attempted}, probes={len(probes)}, skipped: {reasons_skipped}")

    # Write probe FASTA
    probes_fasta = os.path.join(output_dir, 'HCV_probes_consecutive.fasta')
    with open(probes_fasta, 'w') as f:
        for i, (seq_id, probe, pos) in enumerate(probe_records):
            f.write(f">{seq_id}_probe{i+1}_pos{pos}\n{probe}\n")

    # Write summary CSV
    summary_df = pd.DataFrame(summary_rows)
    summary_csv = os.path.join(output_dir, 'probe_design_summary_consecutive.csv')
    summary_df.to_csv(summary_csv, index=False)

    print(f"\nDesigned consecutive probes for {len(summary_rows)} sequences.")
    print(f"Total probes generated: {len(probe_records)}")
    print(f"Probe FASTA: {probes_fasta}")
    print(f"Summary CSV: {summary_csv}")

if __name__ == "__main__":
    import sys
    if len(sys.argv) != 3:
        print("Usage: python3 design_probes_v2.py <input_fasta> <output_dir>")
        sys.exit(1)
    design_probes_consecutive(sys.argv[1], sys.argv[2])
