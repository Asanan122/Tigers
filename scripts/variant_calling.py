#!/usr/bin/env python3

import os
import subprocess
import logging
from pathlib import Path

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

# Define paths
GATK = "gatk-4.4.0.0/gatk"
REFERENCE = "data/reference/panthera_tigris.fasta"
RESULTS_DIR = Path("results")
ALIGNMENT_DIR = RESULTS_DIR / "alignment"
VARIANT_DIR = RESULTS_DIR / "variants"

def run_command(cmd):
    """Run a shell command and log its output."""
    logger.info(f"Running command: {' '.join(cmd)}")
    try:
        subprocess.run(cmd, check=True)
    except subprocess.CalledProcessError as e:
        logger.error(f"Command failed: {e}")
        raise

def mark_duplicates(input_bam, output_bam, metrics_file):
    """Mark duplicates in BAM file."""
    cmd = [
        GATK, "MarkDuplicates",
        "-I", input_bam,
        "-O", output_bam,
        "-M", metrics_file,
        "--CREATE_INDEX", "true"
    ]
    run_command(cmd)

def add_read_groups(input_bam, output_bam, sample_name):
    """Add read groups to BAM file."""
    cmd = [
        GATK, "AddOrReplaceReadGroups",
        "-I", input_bam,
        "-O", output_bam,
        "-RGID", "1",
        "-RGLB", "lib1",
        "-RGPL", "illumina",
        "-RGPU", "unit1",
        "-RGSM", sample_name
    ]
    run_command(cmd)

def haplotype_caller(input_bam, output_gvcf):
    """Run GATK HaplotypeCaller to generate GVCF."""
    cmd = [
        GATK, "HaplotypeCaller",
        "-R", REFERENCE,
        "-I", input_bam,
        "-O", output_gvcf,
        "--emit-ref-confidence", "GVCF",
        "--minimum-mapping-quality", "20",
        "--base-quality-score-threshold", "20"
    ]
    run_command(cmd)

def combine_gvcfs(gvcfs, output_vcf):
    """Combine GVCF files."""
    cmd = [
        GATK, "CombineGVCFs",
        "-R", REFERENCE,
        "-O", output_vcf
    ]
    for gvcf in gvcfs:
        cmd.extend(["-V", str(gvcf)])
    run_command(cmd)

def genotype_gvcfs(input_vcf, output_vcf):
    """Genotype combined GVCF file."""
    cmd = [
        GATK, "GenotypeGVCFs",
        "-R", REFERENCE,
        "-V", input_vcf,
        "-O", output_vcf
    ]
    run_command(cmd)

def filter_variants(input_vcf, output_vcf):
    """Apply variant filtering."""
    cmd = [
        GATK, "VariantFiltration",
        "-R", REFERENCE,
        "-V", input_vcf,
        "-O", output_vcf,
        "--filter-expression", "QD < 2.0 || FS > 60.0 || MQ < 40.0",
        "--filter-name", "basic_filters"
    ]
    run_command(cmd)

def index_bam(input_bam):
    """Index BAM file."""
    cmd = [
        "samtools",
        "index",
        input_bam
    ]
    run_command(cmd)

def main():
    # Create necessary directories
    VARIANT_DIR.mkdir(exist_ok=True)
    
    # Process each sample
    samples = ["caspian", "amur"]
    processed_bams = []
    gvcfs = []
    
    # Check existing files and process only what's needed
    for sample in samples:
        input_bam = ALIGNMENT_DIR / f"{sample}.sorted.bam"
        dedup_bam = VARIANT_DIR / f"{sample}.dedup.bam"
        with_rg_bam = VARIANT_DIR / f"{sample}.with_rg.bam"
        gvcf = VARIANT_DIR / f"{sample}.g.vcf.gz"
        
        # Skip if GVCF already exists
        if gvcf.exists():
            logger.info(f"GVCF for {sample} already exists, skipping processing")
            gvcfs.append(gvcf)
            continue
            
        # Mark duplicates if needed
        if not dedup_bam.exists():
            logger.info(f"Marking duplicates for {sample}")
            mark_duplicates(
                str(input_bam),
                str(dedup_bam),
                str(VARIANT_DIR / f"{sample}.metrics.txt")
            )
        
        # Add read groups if needed
        if not with_rg_bam.exists():
            logger.info(f"Adding read groups for {sample}")
            add_read_groups(
                str(dedup_bam),
                str(with_rg_bam),
                sample
            )
            index_bam(str(with_rg_bam))
        
        # HaplotypeCaller
        logger.info(f"Running HaplotypeCaller for {sample}")
        haplotype_caller(str(with_rg_bam), str(gvcf))
        
        processed_bams.append(with_rg_bam)
        gvcfs.append(gvcf)
    
    # Combine GVCFs
    combined_vcf = VARIANT_DIR / "combined.g.vcf.gz"
    if not combined_vcf.exists():
        logger.info("Combining GVCFs")
        combine_gvcfs([str(g) for g in gvcfs], str(combined_vcf))
    
    # Genotype GVCFs
    genotyped_vcf = VARIANT_DIR / "genotyped.vcf.gz"
    if not genotyped_vcf.exists():
        logger.info("Genotyping GVCFs")
        genotype_gvcfs(str(combined_vcf), str(genotyped_vcf))
    
    # Filter variants
    final_vcf = VARIANT_DIR / "final.vcf.gz"
    if not final_vcf.exists():
        logger.info("Filtering variants")
        filter_variants(str(genotyped_vcf), str(final_vcf))
    
    logger.info("Variant calling pipeline completed successfully!")

if __name__ == "__main__":
    main() 