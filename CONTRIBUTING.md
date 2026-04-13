# Contributing & running from scratch

This guide is for **anyone** (including future you) who needs to run or extend the pipeline. Run commands from the **repository root** unless noted. 

## Prerequisites

- **macOS or Linux** (paths and examples assume Unix shells).
- **Python 3.11+** recommended.
- **Git**, **BWA**, **samtools**, **GATK** (or local `gatk-4.4.0.0/`), **bcftools** (for `call_variants.sh`), **FastQC**, **Java** (SnpEff), and **fasterq-dump** (SRA toolkit) on `PATH` as needed.
- **Disk:** tens of GB free for references, trimmed FASTQ, BAMs, and VCFs.

## 1. Clone and Python env

```bash
git clone https://github.com/Asanan122/Tigers.git
cd Tigers
python3 -m venv venv
source venv/bin/activate
pip install -r requirements.txt
```

## 2. Reference and reads

- **Reference:** use `scripts/download_reference.py` or place `data/reference/panthera_tigris.fasta` (and indexes) per your workflow.
- **Reads:** `scripts/download_genomes.py` pulls SRA data; expect **hours** and large downloads depending on bandwidth.

**Rough runtime:** download + `fasterq-dump` — **hours to overnight** per sample (network and I/O bound).

## 3. Preprocessing

```bash
python scripts/preprocess.py
```

**Outputs:** `results/preprocessing/` (FastQC, trimmed reads).  
**Rough runtime:** **~1–4+ hours** per sample depending on depth and CPU.

## 4. Alignment

```bash
bash scripts/align_reads.sh
# or follow scripts/align_and_call.py
```

**Outputs:** `results/alignment/*.sorted.bam` (+ `.bai`).  
**Rough runtime:** **many hours per sample** for full WGS on a single machine.

## 5. Variant calling (choose a track)

**GATK-style (primary):**

```bash
bash scripts/prepare_resources.sh   # reference dict / faidx if missing
python scripts/variant_calling.py     # or scripts/variant_calling.sh
```

**bcftools (lighter / cross-check):**

```bash
bash scripts/call_variants.sh
```

**Outputs:** VCFs under `data/*/aligned/` (bcftools path) or `results/variants/` (GATK path — follow your script outputs). Merge to a joint VCF (e.g. `results/variants/final.vcf.gz`) using your chosen workflow.

**Rough runtime:** **hours to days** for BQSR + HaplotypeCaller-style work at WGS depth; bcftools track is usually faster but less feature-rich.

## 6. Analysis and figures

```bash
python scripts/analyze_and_visualize.py
python scripts/detailed_analysis.py
```

**Outputs:** `results/analysis/*.png`, `variant_statistics.txt`, `detailed_statistics.txt`, large CSVs.

**Rough runtime:** **~tens of minutes to a few hours** depending on VCF size and machine (I/O + Python).

## 7. SnpEff (optional)

```bash
python scripts/setup_snpeff_db.py
python scripts/annotate_variants.py
```

Ensure VCF contig names match the GFF/reference used for SnpEff.

## 8. Pairwise / gene-level scripts

```bash
python scripts/compare_genomes.py
python scripts/compare_tiger_genomes.py
```

`compare_tiger_genomes.py` writes **`variant_analysis.log`** in the **current working directory** if you run it from elsewhere — prefer repo root.

## Pull requests & issues

- **Small, focused changes** (one concern per PR) are easiest to review.
- **Document** new scripts in `README.md` and, if they change outputs, update **`RESULTS.md`**.
- Prefer **pinned** dependencies if you add Python packages (`pip freeze` after intentional upgrades).

## License

Add a `LICENSE` file if you intend open redistribution; until then, assume **all rights reserved** unless you state otherwise in the repository.
