# Results summary — Caspian vs Amur nuclear WGS

This file summarizes **headline outputs** from completed pipeline runs so the repository remains interpretable without cloning multi-gigabyte `data/` or `results/`. Raw VCFs, BAMs, and large CSVs stay local (see `.gitignore`). 

## Study design (brief)

- **Samples:** one Caspian (SRR18572400) and one Amur (SRR31485304) whole-genome resequencing dataset, aligned to a *Panthera tigris* reference.
- **Scope:** nuclear genome variation **between two individuals** (not a full population study). Statistics below describe **joint VCF** sites and genotype presence patterns.

## Headline variant statistics

**Primary quantitative report:** `results/analysis_summary.txt` — regenerate anytime with:

```bash
python scripts/analysis_summary_report.py
```

That file merges a full pass over the joint VCF (QUAL, contigs, diploid genotype patterns) with cached `results/analysis/*.txt` outputs where present.

Legacy table (from `scripts/analyze_and_visualize.py`; sample column labels may not match VCF order — see note inside `analysis_summary.txt`):

| Metric | Value | Notes |
|--------|------:|-------|
| Total variant rows (joint VCF) | 2,169,502 | Rows in merged call set |
| Sites with a call in **both** samples | 73,119 | Both genotypes non-missing at the site |
| Sites called **Caspian only** (Amur missing) | 2,008,535 | Often reflects coverage / QC asymmetry, not only “private subspecies SNPs” |
| Sites called **Amur only** (Caspian missing) | 87,848 | Same caveat |
| Opposite homozygous genotypes (0 vs 2) | 194 | Strongest simple pairwise contrast among called genotypes |

### Variant class breakdown

From `results/analysis/detailed_statistics.txt` (pandas summary on the joint table):

| Class | Approx. count |
|-------|----------------:|
| SNPs | 1,919,574 |
| Insertions | 114,467 |
| Deletions | 135,398 |

## Chromosomes and density

Per-contig counts and density plots are generated under **`results/analysis/`** when you run `detailed_analysis.py` / `analyze_and_visualize.py`, e.g.:

- `chromosome_distribution.png` — variant counts by contig  
- `variant_density.png` — density along positions (joint VCF)

Re-run those scripts after regenerating `final.vcf.gz` to refresh figures.

## Functional annotation (SnpEff)

In the exported gene summaries from earlier runs, **most records lacked resolved gene names in the ANN field** (`Unknown` in `results/analysis/gene_summary.csv`), and **high-impact** rows were essentially empty until annotation and contig naming are fully harmonized with the reference GFF. **Do not over-interpret functional categories** until SnpEff is run end-to-end with a DB that matches VCF contig names.

## Narrative write-up

A longer draft discussion lives in **`results/manuscript.md`** (tracked in git when present). Treat interpretation there as **hypothesis-generating** until sample sizes and annotation are strengthened.

## Reproducing this document’s numbers

1. Produce or obtain `results/variants/final.vcf.gz` (joint, normalized to your convention).  
2. From repo root: `python scripts/analyze_and_visualize.py`  
3. Read `results/analysis/variant_statistics.txt` and `detailed_statistics.txt`.
