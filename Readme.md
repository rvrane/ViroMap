# ViroMap

> **Targeted virome capture panel design pipeline** for hybrid-capture NGS enrichment of vertebrate DNA and RNA viruses from clinical and environmental samples.

![Language](https://img.shields.io/badge/language-Python%20%7C%20Shell-blue)
![Platform](https://img.shields.io/badge/platform-Linux%20%7C%20HPC-lightgrey)
![License](https://img.shields.io/badge/license-MIT-green)

---

## Overview

ViroMap is a computational probe design pipeline for constructing targeted virome capture panels. The panel targets vertebrate-infecting DNA and RNA viruses with a non-redundant probe set suitable for solution-based hybrid capture workflows. When coupled with standard metagenomic sequencing, this design enables several-hundred-fold enrichment of viral read fractions and markedly improved genome breadth and depth compared to unenriched libraries.

**Key applications:**
- Comprehensive virome profiling in clinical, surveillance, and research settings
- Sensitive detection of low-abundance or mixed viral infections in high background host material
- Characterization of both known viruses and divergent variants with substantial sequence mismatch to the probe design

---

## Repository Structure

```
ViroMap/
├── cluster_sequences_v2.sh         # CD-HIT based viral sequence clustering (v2)
├── cluster_sequences_v3.sh         # Improved clustering with updated thresholds (v3)
├── design_probes_v3.py             # Probe tiling and design logic from clustered sequences
├── filter_probes_v1.sh             # Post-design probe filtering (GC content, complexity)
├── remove_degenerate_bases.py      # Removes probes containing degenerate IUPAC bases
├── validate_probes_blast_v2.sh     # BLAST-based off-target and specificity validation
└── Degenerate_removal_Summary.md   # Summary report of degenerate base filtering step
```

---

## Pipeline Workflow

```
Viral Reference Sequences (NCBI/ViPR)
            │
            ▼
  1. Cluster Sequences         (cluster_sequences_v3.sh)
     └─ CD-HIT-EST clustering to reduce redundancy
            │
            ▼
  2. Design Probes             (design_probes_v3.py)
     └─ Sliding-window tiling of representative sequences
            │
            ▼
  3. Filter Probes             (filter_probes_v1.sh)
     └─ GC-content, low-complexity, and length filters
            │
            ▼
  4. Remove Degenerate Bases   (remove_degenerate_bases.py)
     └─ Strict IUPAC base validation
            │
            ▼
  5. Validate Specificity      (validate_probes_blast_v2.sh)
     └─ BLAST against host genome to remove off-target probes
            │
            ▼
       Final Probe Panel (.fasta)
```

---

## Requirements

| Tool | Version | Purpose |
|------|---------|-------- |
| Python | ≥ 3.8 | Probe design and degenerate base removal |
| CD-HIT-EST | ≥ 4.8 | Sequence clustering |
| BLAST+ | ≥ 2.11 | Off-target probe validation |
| Bash | ≥ 4.0 | Shell scripts |

---

## Installation

```bash
# Clone the repository
git clone https://github.com/rvrane/ViroMap.git
cd ViroMap

# Make shell scripts executable
chmod +x *.sh

# Install Python dependencies
pip install biopython
```

---

## Usage

### Step 1: Cluster viral reference sequences
```bash
bash cluster_sequences_v3.sh \
  --input viral_references.fasta \
  --output clustered_representatives.fasta \
  --identity 0.90
```

### Step 2: Design probes
```bash
python design_probes_v3.py \
  --input clustered_representatives.fasta \
  --probe_length 120 \
  --tiling_offset 60 \
  --output probes_raw.fasta
```

### Step 3: Filter and clean probes
```bash
bash filter_probes_v1.sh --input probes_raw.fasta --output probes_filtered.fasta
python remove_degenerate_bases.py --input probes_filtered.fasta --output probes_clean.fasta
```

### Step 4: Validate specificity
```bash
bash validate_probes_blast_v2.sh \
  --input probes_clean.fasta \
  --host_db /path/to/hg38_blast_db \
  --output probes_validated.fasta
```

---

## Output

| File | Description |
|------|-------------|
| `clustered_representatives.fasta` | Non-redundant viral sequence representatives |
| `probes_raw.fasta` | All tiled probes before filtering |
| `probes_filtered.fasta` | GC and complexity filtered probes |
| `probes_clean.fasta` | Degenerate-base-free probe set |
| `probes_validated.fasta` | Final panel after host off-target removal |
| `Degenerate_removal_Summary.md` | QC report for degenerate base removal |

---

## Notes & Limitations

- Probe sensitivity depends on the completeness of the input viral reference database
- CD-HIT clustering identity threshold (`--identity`) should be tuned based on target viral diversity
- Host off-target validation requires a pre-built BLAST database for the relevant host genome
- Designed for Linux/HPC environments; not tested on Windows natively

---

## Future Work

- [ ] Integration with Snakemake or Nextflow for end-to-end workflow management
- [ ] Docker/Conda containerization
- [ ] Automated NCBI viral sequence fetching module
- [ ] Support for RNA virus-specific uracil handling

---

## Author

**Rugved Rane**  
Senior Bioinformatician | Clinical Genomics  
[GitHub](https://github.com/rvrane)

---

## License

MIT License — see [LICENSE](LICENSE) for details.
