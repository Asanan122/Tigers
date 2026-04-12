# Tiger genome analysis (Caspian vs Amur)

Whole-genome resequencing workflow comparing **Caspian** (*Panthera tigris virgata*, SRR18572400) and **Amur** (*Panthera tigris altaica*, SRR31485304) to a *Panthera tigris* reference: QC and trimming, alignment, variant calling, annotation, and downstream summaries and figures.

## Layout

```
.
├── data/
│   ├── amur/                 # Amur FASTQ + per-sample aligned/ outputs
│   ├── caspian/              # Caspian FASTQ + per-sample aligned/ outputs
│   └── reference/            # Reference FASTA, indexes, GFF, NCBI dataset payload, snpeff_db/
├── gatk-4.4.0.0/             # Bundled GATK (used by variant_calling / prepare_resources)
├── tools/
│   └── snpEff/               # Local only (gitignored): install via `annotate_variants.py` or SnpEff docs
├── results/
│   ├── preprocessing/        # FastQC, trimmed reads
│   ├── alignment/            # Sorted BAMs + .bai
│   ├── variants/             # VCFs and related outputs
│   ├── analysis/             # CSV summaries, plots, statistics
│   ├── comparison/plots/     # Reserved for comparison figures
│   ├── figures/
│   ├── manuscript.md         # Draft write-up
│   └── *.csv / *.png         # Some top-level comparison outputs
├── scripts/                  # All automation (see below)
├── requirements.txt
└── README.md
```

## Python environment

Create a virtual environment and install dependencies (the repo does **not** ship `venv/`):

```bash
cd /path/to/Tigers
python3 -m venv venv
source venv/bin/activate   # Windows: venv\Scripts\activate
pip install -r requirements.txt
```

## External tools

Install on your `PATH` as needed:

| Tool | Typical use |
|------|-------------|
| **BWA** | Read alignment |
| **samtools** | BAM/CRAM, indexing, faidx |
| **GATK** | Bundled at `gatk-4.4.0.0/gatk` (used by `variant_calling.py` / shell scripts); or install elsewhere and edit script paths |
| **bcftools** | `call_variants.sh` (mpileup + call) |
| **FastQC** | QC reports |
| **Java** | SnpEff (`java -jar tools/snpEff/snpEff.jar`) |
| **fasterq-dump** | `download_genomes.py` (SRA) |

## Scripts (by role)

| Script | Role |
|--------|------|
| `download_reference.py`, `download_genomes.py` | Fetch reference / SRA data |
| `preprocess.py` | QC and trimming |
| `align_reads.sh`, `align_and_call.py` | Alignment |
| `variant_calling.py`, `variant_calling.sh` | GATK-oriented variant calling |
| `call_variants.sh` | bcftools mpileup/call per sample |
| `prepare_resources.sh` | Reference dict / index via GATK + samtools |
| `compare_genomes.py`, `compare_tiger_genomes.py` | Compare calls; gene-level stats + plots (writes `variant_analysis.log` in CWD when run) |
| `analyze_and_visualize.py`, `detailed_analysis.py`, `analyze_variants.sh` | Summaries and plots under `results/analysis/` |
| `annotate_variants.py` | Download/configure SnpEff and annotate VCFs |
| `setup_snpeff_db.py` | Copy reference files into `tools/snpEff/data/panthera_tigris/` and write `snpEff.config` |
| `convert_gff.py`, `convert_gff_chromosomes.py` | GFF helpers for SnpEff / reference |

Run Python scripts from the **repository root** so relative paths resolve.

## Pipeline status (high level)

- Preprocessing (FastQC, trimming) and alignment outputs are present under `results/preprocessing/` and `results/alignment/`.
- Variant calling and downstream analysis artifacts live under `results/variants/` and `results/analysis/`; see `results/manuscript.md` for a narrative summary of earlier results.

## Cleanup notes

The following are intentionally **not** kept in the repo (regenerate or reinstall as needed):

- **`venv/`** — recreate with `python3 -m venv venv` and `pip install -r requirements.txt`.
- **`fasterq.tmp.*/`** — temporary SRA/fasterq-dump scratch; safe to delete if it reappears.
- **`tools/snpEff.zip`** — redundant once `tools/snpEff/` is unpacked; `annotate_variants.py` can re-fetch.
- **SnpEff `examples/` and `galaxy/`** — removed to save space; not required to run annotation on your VCFs.
- **`tools/clinEff/`** — unused by this pipeline.

Large directories stay local only (see `.gitignore`): `data/`, most of `results/`, bundled GATK (`gatk-4.4.0.0/` and `gatk-4.4.0.0.zip`), and all of `tools/snpEff/` (reinstall after clone). Only `results/manuscript.md` and `results/analysis_summary.txt` are tracked for GitHub.

### GitHub

After clone, restore reference reads, alignments, and outputs locally (re-run download / pipeline scripts), then run `setup_snpeff_db.py` and install GATK/SnpEff as in this README.

## Requirements

See `requirements.txt` (numpy, pandas, matplotlib, seaborn, cyvcf2, pysam, gffutils, requests, tqdm).
