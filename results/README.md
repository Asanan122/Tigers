# `results/` — what is in Git vs local-only

Most pipeline outputs (BAMs, VCFs, large CSVs, plots) are **gitignored** so clones stay small. Regenerate them with the root **[CONTRIBUTING.md](../CONTRIBUTING.md)** workflow.

## Files tracked in this repository

| File | Description |
|------|-------------|
| **[manuscript.md](manuscript.md)** | Draft narrative of variant findings and interpretation |
| **[analysis_summary.txt](analysis_summary.txt)** | Quantitative joint-VCF summary; refresh with `python scripts/analysis_summary_report.py` from repo root |

On GitHub, open these paths directly if the folder view is sparse: they are normal tracked files under `results/`.
