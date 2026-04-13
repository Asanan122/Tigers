"""
Build results/analysis_summary.txt from the joint VCF plus optional pipeline artifacts.
Run from repository root: python scripts/analysis_summary_report.py
"""

from __future__ import annotations

import os
import sys
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path
from numbers import Real
from typing import Any, Dict, List, Optional, Tuple

import numpy as np

try:
    import cyvcf2
except ImportError as e:
    raise SystemExit(
        "cyvcf2 is not installed for this Python interpreter.\n"
        "  Use the project venv (recommended):  source venv/bin/activate && pip install -r requirements.txt\n"
        "  Then run:  python scripts/analysis_summary_report.py\n"
        "Avoid `pip install` into base Anaconda — it can upgrade numpy and break scipy/numba/gensim."
    ) from e


def _root() -> Path:
    return Path(__file__).resolve().parent.parent


def _gt_alleles(g) -> Tuple[int, int]:
    return int(g[0]), int(g[1])


def _is_called(g) -> bool:
    a, b = _gt_alleles(g)
    return a >= 0 and b >= 0


def _is_hom_ref(g) -> bool:
    a, b = _gt_alleles(g)
    return a == 0 and b == 0


def _is_hom_alt(g) -> bool:
    """Homozygous non-reference (biallelic-first-alt assumption)."""
    a, b = _gt_alleles(g)
    if a < 0 or b < 0:
        return False
    return a == b and a != 0


def _is_het(g) -> bool:
    return _is_called(g) and not _is_hom_ref(g) and not _is_hom_alt(g)


def _variant_class(ref: str, alt0: Optional[str]) -> str:
    if not alt0:
        return "other"
    if len(ref) == 1 and len(alt0) == 1:
        return "SNP"
    if len(ref) > len(alt0):
        return "deletion"
    if len(ref) < len(alt0):
        return "insertion"
    return "complex"


def _passes_filter(rec) -> bool:
    flt = rec.FILTER
    if not flt:
        return True
    if isinstance(flt, str):
        return flt in ("PASS", ".", "")
    try:
        return list(flt) == [] or "PASS" in flt
    except TypeError:
        return True


def compute_joint_vcf_stats(vcf_path: Path) -> Dict[str, Any]:
    """Single pass over joint VCF; returns serializable summary dict."""
    vcf = cyvcf2.VCF(str(vcf_path))
    samples: List[str] = list(vcf.samples)
    if len(samples) < 2:
        vcf.close()
        raise ValueError(f"Expected ≥2 samples in {vcf_path}, found {samples}")

    n = 0
    by_chrom: Dict[str, int] = defaultdict(int)
    snp = indel_del = ins = complex_ = 0
    pass_n = 0
    quals: List[float] = []

    both_called = 0
    s0_only = s1_only = 0
    neither_called = 0
    hom_ref_both = hom_alt_both = 0
    opposite_hom = 0
    same_genotype = 0
    both_het = 0

    # Per-sample called counts (any non-missing genotype)
    called = [0, 0]
    hom_ref = [0, 0]
    hom_alt = [0, 0]
    het = [0, 0]

    idx0, idx1 = 0, 1
    name0, name1 = samples[idx0], samples[idx1]

    for rec in vcf:
        n += 1
        by_chrom[rec.CHROM] += 1
        if _passes_filter(rec):
            pass_n += 1

        q = rec.QUAL
        if q is not None and not (isinstance(q, float) and np.isnan(q)):
            quals.append(float(q))

        alt0 = rec.ALT[0] if rec.ALT else None
        alt_s = alt0.decode() if isinstance(alt0, (bytes, bytearray)) else (alt0 or "")
        ref_s = rec.REF.decode() if isinstance(rec.REF, (bytes, bytearray)) else rec.REF
        vc = _variant_class(ref_s, alt_s)
        if vc == "SNP":
            snp += 1
        elif vc == "deletion":
            indel_del += 1
        elif vc == "insertion":
            ins += 1
        else:
            complex_ += 1

        g0 = rec.genotypes[idx0]
        g1 = rec.genotypes[idx1]

        for i, g in enumerate((g0, g1)):
            if _is_called(g):
                called[i] += 1
                if _is_hom_ref(g):
                    hom_ref[i] += 1
                elif _is_hom_alt(g):
                    hom_alt[i] += 1
                elif _is_het(g):
                    het[i] += 1

        c0, c1 = _is_called(g0), _is_called(g1)
        if c0 and c1:
            both_called += 1
            if _is_hom_ref(g0) and _is_hom_ref(g1):
                hom_ref_both += 1
            elif _is_hom_alt(g0) and _is_hom_alt(g1):
                hom_alt_both += 1
            elif _is_het(g0) and _is_het(g1):
                both_het += 1
            ga = sorted(_gt_alleles(g0))
            gb = sorted(_gt_alleles(g1))
            if ga == gb:
                same_genotype += 1
            if (_is_hom_ref(g0) and _is_hom_alt(g1)) or (_is_hom_alt(g0) and _is_hom_ref(g1)):
                opposite_hom += 1
        elif c0 and not c1:
            s0_only += 1
        elif c1 and not c0:
            s1_only += 1
        else:
            neither_called += 1

    vcf.close()

    qarr = np.array(quals, dtype=np.float64) if quals else np.array([], dtype=np.float64)

    top_chrom = sorted(by_chrom.items(), key=lambda x: -x[1])[:15]

    return {
        "vcf_path": str(vcf_path),
        "sample_0": name0,
        "sample_1": name1,
        "total_records": n,
        "n_contigs_with_variants": len(by_chrom),
        "snp": snp,
        "deletion": indel_del,
        "insertion": ins,
        "complex": complex_,
        "pass_records": pass_n,
        "pass_fraction": pass_n / n if n else 0.0,
        "qual_mean": float(qarr.mean()) if len(qarr) else float("nan"),
        "qual_median": float(np.median(qarr)) if len(qarr) else float("nan"),
        "qual_std": float(qarr.std()) if len(qarr) else float("nan"),
        "qual_p10": float(np.percentile(qarr, 10)) if len(qarr) else float("nan"),
        "qual_p90": float(np.percentile(qarr, 90)) if len(qarr) else float("nan"),
        "both_called": both_called,
        "only_sample0": s0_only,
        "only_sample1": s1_only,
        "neither_called": neither_called,
        "hom_ref_both": hom_ref_both,
        "hom_alt_both": hom_alt_both,
        "both_het": both_het,
        "same_genotype_both_called": same_genotype,
        "opposite_homozygous": opposite_hom,
        "called_counts": {"sample0": called[0], "sample1": called[1]},
        "hom_ref_counts": {"sample0": hom_ref[0], "sample1": hom_ref[1]},
        "hom_alt_counts": {"sample0": hom_alt[0], "sample1": hom_alt[1]},
        "het_counts": {"sample0": het[0], "sample1": het[1]},
        "top_contigs": top_chrom,
    }


def _read_text_if_exists(path: Path) -> Optional[str]:
    if not path.is_file():
        return None
    return path.read_text(encoding="utf-8", errors="replace").strip()


def write_comprehensive_analysis_summary(
    joint_vcf: Optional[Path] = None,
    gene_summary: Optional[Dict[str, Any]] = None,
    output_path: Optional[Path] = None,
) -> Path:
    root = _root()
    os.chdir(root)
    joint_vcf = joint_vcf or root / "results" / "variants" / "final.vcf.gz"
    output_path = output_path or root / "results" / "analysis_summary.txt"

    lines: List[str] = []
    ts = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC")
    lines.append("=" * 78)
    lines.append("TIGER WGS — QUANTITATIVE ANALYSIS SUMMARY")
    lines.append("=" * 78)
    lines.append(f"Generated: {ts}")
    lines.append("")

    if not joint_vcf.is_file():
        lines.append(f"[WARN] Joint VCF not found: {joint_vcf}")
        lines.append("       Place final.vcf.gz or pass joint_vcf=... to regenerate stats.")
        lines.append("")
    else:
        lines.append("1. JOINT CALL SET (single pass over VCF)")
        lines.append("-" * 78)
        st = compute_joint_vcf_stats(joint_vcf)
        lines.append(f"  VCF path: {st['vcf_path']}")
        lines.append(f"  Sample index 0: {st['sample_0']}")
        lines.append(f"  Sample index 1: {st['sample_1']}")
        lines.append("")
        lines.append(f"  Total variant records: {st['total_records']:,}")
        lines.append(f"  Contigs with ≥1 record: {st['n_contigs_with_variants']:,}")
        lines.append(f"  SNPs (len(ref)==len(alt)==1): {st['snp']:,}")
        lines.append(f"  Deletions: {st['deletion']:,}")
        lines.append(f"  Insertions: {st['insertion']:,}")
        lines.append(f"  Other / complex length classes: {st['complex']:,}")
        lines.append("")
        lines.append("2. QUALITY AND FILTER")
        lines.append("-" * 78)
        lines.append(f"  Records with PASS (or empty FILTER): {st['pass_records']:,} ({100*st['pass_fraction']:.2f}%)")
        lines.append(
            f"  QUAL (numeric): mean={st['qual_mean']:.2f}, median={st['qual_median']:.2f}, "
            f"stdev={st['qual_std']:.2f}"
        )
        lines.append(
            f"  QUAL percentiles: p10={st['qual_p10']:.2f}, p90={st['qual_p90']:.2f}"
        )
        lines.append("")
        lines.append("3. GENOTYPE PATTERNS (diploid; labels follow VCF sample order)")
        lines.append("-" * 78)
        lines.append(
            f"  Sites with a call in BOTH samples: {st['both_called']:,} "
            f"({100*st['both_called']/st['total_records']:.2f}% of rows)"
        )
        lines.append(
            f"  Called only in {st['sample_0']}: {st['only_sample0']:,} "
            f"({100*st['only_sample0']/st['total_records']:.2f}%)"
        )
        lines.append(
            f"  Called only in {st['sample_1']}: {st['only_sample1']:,} "
            f"({100*st['only_sample1']/st['total_records']:.2f}%)"
        )
        lines.append(
            f"  Missing / uncalled in BOTH: {st['neither_called']:,} "
            f"({100*st['neither_called']/st['total_records']:.2f}%)"
        )
        lines.append("")
        lines.append("  Among sites with BOTH samples called:")
        lines.append(f"    Same genotype (allele pair): {st['same_genotype_both_called']:,}")
        lines.append(f"    Both homozygous reference: {st['hom_ref_both']:,}")
        lines.append(f"    Both homozygous alternate: {st['hom_alt_both']:,}")
        lines.append(f"    Both heterozygous: {st['both_het']:,}")
        lines.append(
            f"    Opposite homozygous (ref vs alt), diploid-aware: {st['opposite_homozygous']:,} "
            f"— contrast where one sample is hom-REF and the other hom-ALT (both called)"
        )
        lines.append(
            "    (Legacy `analyze_and_visualize.py` reported 'Homozygous Differences: 194' using only "
            "the first cyvcf2 allele index per sample and an assumed sample order; use the line above "
            "for diploid biology.)"
        )
        lines.append("")
        lines.append("  Per-sample genotype counts (rows where that sample is called):")
        lines.append(
            f"    {st['sample_0']}: called={st['called_counts']['sample0']:,}, "
            f"hom_ref={st['hom_ref_counts']['sample0']:,}, "
            f"hom_alt={st['hom_alt_counts']['sample0']:,}, "
            f"het={st['het_counts']['sample0']:,}"
        )
        lines.append(
            f"    {st['sample_1']}: called={st['called_counts']['sample1']:,}, "
            f"hom_ref={st['hom_ref_counts']['sample1']:,}, "
            f"hom_alt={st['hom_alt_counts']['sample1']:,}, "
            f"het={st['het_counts']['sample1']:,}"
        )
        lines.append("")
        lines.append("4. TOP CONTIGS BY VARIANT COUNT")
        lines.append("-" * 78)
        for chrom, cnt in st["top_contigs"]:
            lines.append(f"  {chrom}: {cnt:,}")
        lines.append("")
        lines.append("OBSERVATION — asymmetry in 'called in one sample only':")
        lines.append(
            "  Large one-sided missingness usually reflects coverage, batch, or joint"
        )
        lines.append(
            "  genotyping rather than biological 'private subspecies SNPs' alone."
        )
        lines.append(
            "  Interpret opposite-homozygous counts as pairwise contrast among called sites."
        )
        lines.append("")

    aux_vs = root / "results" / "analysis" / "variant_statistics.txt"
    aux_ds = root / "results" / "analysis" / "detailed_statistics.txt"
    txt_vs = _read_text_if_exists(aux_vs)
    txt_ds = _read_text_if_exists(aux_ds)
    if txt_vs or txt_ds:
        lines.append("5. CACHED TEXT SUMMARIES (from prior analyze_and_visualize / detailed_analysis)")
        lines.append("-" * 78)
        lines.append(
            "  NOTE: `analyze_and_visualize.py` hard-codes genotype column 0 = 'CASPIAN' and 1 = 'AMUR'."
        )
        lines.append(
            "  This project's `final.vcf.gz` sample order is ['amur','caspian'], so pasted labels in"
        )
        lines.append(
            "  variant_statistics.txt do NOT match VCF order — trust Section 3 for named-sample counts."
        )
        lines.append("")
        if txt_vs:
            lines.append(aux_vs.name + ":")
            lines.append(txt_vs)
            lines.append("")
        if txt_ds:
            lines.append(aux_ds.name + ":")
            lines.append(txt_ds)
            lines.append("")

    lines.append("6. GENE-LEVEL OVERLAP (optional; separate per-sample VCF + GFF path)")
    lines.append("-" * 78)
    if gene_summary:
        for key, value in sorted(gene_summary.items()):
            if isinstance(value, Real) and not isinstance(value, bool):
                lines.append(f"  {key.replace('_', ' ').title()}: {float(value):.4f}")
            else:
                lines.append(f"  {key.replace('_', ' ').title()}: {value}")
    else:
        lines.append(
            "  Not merged into this run. For gene windows vs per-sample VCFs, run:"
        )
        lines.append("    python scripts/compare_tiger_genomes.py")
        lines.append(
            "  (Requires GFF/VCF contig naming consistency or chromosome_mapping.txt fixes.)"
        )
    lines.append("")
    lines.append("=" * 78)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return output_path


def main() -> None:
    root = _root()
    os.chdir(root)
    out = write_comprehensive_analysis_summary()
    print(f"Wrote {out}")


if __name__ == "__main__":
    main()
