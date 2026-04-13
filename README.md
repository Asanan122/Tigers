# Tiger conservation genomics — Caspian vs Amur nuclear WGS

## Why this matters

*Panthera tigris* subspecies are among the clearest examples of **intraspecific genomic structure in endangered megafauna**: geographically separated lineages that share a species-level identity but carry distinct allele frequencies, histories of drift and selection, and conservation-relevant diversity. **Understanding what distinguishes living subspecies at the DNA level** is the same class of problem as comparing an extant relative to a **de-extinction or genetic-rescue proxy**: you need reproducible alignment, variant calling, and interpretation against a reference, with explicit uncertainty.

This repository is a **practical whole-genome resequencing workflow** for two public SRA samples, **Caspian** (*P. t. virgata*, SRR18572400) and **Amur** (*P. t. altaica*, SRR31485304), aligned to a *Panthera tigris* reference. It does not claim population-level inference from *N* = 1 per lineage; it **does** provide a defensible pipeline and documented outputs for comparative nuclear variation. See **[RESULTS.md](RESULTS.md)** for headline statistics from completed runs.

## Variant calling: GATK and bcftools

Two callers appear on purpose:

- **`variant_calling.py` / `variant_calling.sh`** — **GATK 4** best-practices style workflow (e.g. MarkDuplicates, BQSR-oriented steps as implemented) for the primary analysis track.
- **`call_variants.sh`** — **bcftools** `mpileup` + `call` for a **lighter-weight, orthogonal** genotype path. That supports **sanity checks**, teaching, or environments where a full GATK run is impractical; discrepancies between engines are expected and worth inspecting at known truth sites or subsampled regions.

Joint downstream summaries (e.g. `analyze_and_visualize.py`) assume a **joint VCF** (e.g. `results/variants/final.vcf.gz`) produced by your chosen path — document which caller generated the file you cite.

## Quick links

| Doc | Purpose |
|-----|---------|
| **[RESULTS.md](RESULTS.md)** | Headline findings and where plots live (visible without large `data/` or `results/` on disk) |
| **[CONTRIBUTING.md](CONTRIBUTING.md)** | End-to-end run order, commands, expected artifacts, rough runtimes |

## Public data sources (recreate from scratch)

All inputs are public; nothing in this workflow depends on private files.

### Whole-genome sequencing reads (NCBI SRA)

| Subspecies (lineage) | Run accession | NCBI SRA record |
|----------------------|---------------|-----------------|
| Caspian (*P. t. virgata*) | **SRR18572400** | [SRR18572400](https://www.ncbi.nlm.nih.gov/sra/SRR18572400) |
| Amur (*P. t. altaica*) | **SRR31485304** | [SRR31485304](https://www.ncbi.nlm.nih.gov/sra/SRR31485304) |

**Fetch locally:** install [SRA Toolkit](https://github.com/ncbi/sra-tools/wiki/01.-Downloading-SRA-Toolkit), then either run `python scripts/download_genomes.py` (uses `fasterq-dump` into `data/caspian/` and `data/amur/`) or manually:

```bash
prefetch SRR18572400 SRR31485304
fasterq-dump --progress --outdir data/caspian SRR18572400
fasterq-dump --progress --outdir data/amur    SRR31485304
```

Paired-end layout and metadata for library prep are on the SRA run pages above.

### Reference genome and index (Ensembl — *Panthera tigris altaica*, PanTig1.0)

The bundled downloader uses **Ensembl release 110** toplevel DNA for **Amur tiger** assembly **PanTig1.0** (same species reference used for read alignment in this project):

- **FASTA (gzip):** [Panthera_tigris_altaica.PanTig1.0.dna.toplevel.fa.gz](http://ftp.ensembl.org/pub/release-110/fasta/panthera_tigris_altaica/dna/Panthera_tigris_altaica.PanTig1.0.dna.toplevel.fa.gz)
- **Species / assembly hub:** [Ensembl *Panthera tigris altaica*](https://www.ensembl.org/Panthera_tigris_altaica/Info/Index)

**Fetch locally:** `python scripts/download_reference.py` writes `data/reference/panthera_tigris.fasta` (decompressed). Then index with BWA / `samtools faidx` / GATK `CreateSequenceDictionary` as in **CONTRIBUTING.md**.

### Annotation (GFF / genes)

Gene models and GFF used in this repo were aligned to the **same reference build** as the FASTA above (check contig names against your VCF before SnpEff). If you need a fresh annotation package, use **Ensembl** for PanTig1.0 or **[NCBI Datasets](https://www.ncbi.nlm.nih.gov/datasets/)** for the matching *Panthera tigris* assembly — **contig naming must match** between FASTA, BAM, VCF, and GFF or downstream scripts will show empty gene overlap.

## Repository layout

```
.
├── data/                    # Local only (gitignored): reads, reference, alignments
├── gatk-4.4.0.0/            # Local GATK bundle (gitignored) or install system-wide
├── tools/snpEff/            # Local SnpEff + DB (gitignored); bootstrap via scripts
├── results/                 # Mostly gitignored; small docs may be tracked
├── scripts/                 # Pipeline entry points
├── requirements.txt         # Fully pinned (pip freeze)
├── RESULTS.md
├── CONTRIBUTING.md
└── README.md
```

## Python environment

```bash
cd /path/to/Tigers
python3 -m venv venv
source venv/bin/activate    # Windows: venv\Scripts\activate
pip install -r requirements.txt
```

**Use this venv for every project command** (`python scripts/...`). If you `pip install` into **base Conda** instead, you can upgrade **numpy** and break other Conda packages (scipy, numba, gensim, etc.). If that already happened, reinstall those packages in Conda or create a fresh Conda env for non-Tigers work.

Dependencies are **pinned to exact versions** for reproducibility. To refresh after upgrading Python or packages intentionally:

```bash
pip install -r requirements.txt
pip freeze > requirements.txt   # only when you mean to change pins project-wide
```

## External tools

| Tool | Role |
|------|------|
| **BWA** | Alignment |
| **samtools** | BAM / CRAM, `faidx` |
| **GATK** | Primary variant pipeline (`gatk-4.4.0.0/gatk` or your install) |
| **bcftools** | `call_variants.sh` |
| **FastQC** | QC reports |
| **Java** | SnpEff |
| **fasterq-dump** | `download_genomes.py` (SRA) |

## Scripts (by role)

| Script | Role |
|--------|------|
| `download_reference.py`, `download_genomes.py` | Reference / SRA fetch |
| `preprocess.py` | QC and trimming |
| `align_reads.sh`, `align_and_call.py` | Alignment |
| `variant_calling.py`, `variant_calling.sh` | GATK-oriented calling |
| `call_variants.sh` | bcftools mpileup / call |
| `prepare_resources.sh` | Sequence dictionary, reference index |
| `compare_genomes.py`, `compare_tiger_genomes.py` | Pairwise comparison, gene-level plots (see `variant_analysis.log` in CWD) |
| `analysis_summary_report.py` | Writes `results/analysis_summary.txt` (joint VCF stats, QUAL, contigs, diploid genotypes) |
| `analyze_and_visualize.py`, `detailed_analysis.py` | Summaries and plots → `results/analysis/` |
| `annotate_variants.py` | SnpEff download / annotate |
| `setup_snpeff_db.py` | Panthera tigris SnpEff DB layout + config |
| `convert_gff.py`, `convert_gff_chromosomes.py` | GFF utilities |

Run Python from the **repository root** unless a script documents otherwise.

## What is (and is not) on GitHub

Large paths are **gitignored** (`data/`, most of `results/`, bundled GATK, `tools/snpEff/`). **RESULTS.md** carries summary statistics so the repository still reads as a **study**, not only a script dump. After clone, restore data locally and re-run or copy outputs per **CONTRIBUTING.md**.


---

*This project is independent work for research and learning; it is not affiliated with Colossal Biosciences. The framing above describes how comparative genomics for endangered felids connects to broader conservation and de-extinction-adjacent science.*
