import os
import subprocess
from pathlib import Path

def run_bwa_index(reference):
    """Index reference genome with BWA"""
    print(f"Indexing reference genome {reference}...")
    cmd = f"bwa index {reference}"
    subprocess.run(cmd, shell=True, check=True)
    print("Reference genome indexing completed.")

def run_bwa_mem(reference, reads1, reads2, output):
    """Align paired-end reads to reference using BWA-MEM"""
    print(f"Aligning reads to reference genome...")
    cmd = f"bwa mem {reference} {reads1} {reads2} | samtools view -bS - | samtools sort -o {output}"
    subprocess.run(cmd, shell=True, check=True)
    print("Alignment completed.")

def run_samtools_index(bam_file):
    """Index BAM file"""
    print(f"Indexing BAM file {bam_file}...")
    cmd = f"samtools index {bam_file}"
    subprocess.run(cmd, shell=True, check=True)
    print("BAM indexing completed.")

def run_bcftools_call(bam_file, output_vcf, reference):
    """Call variants using bcftools"""
    print(f"Calling variants for {bam_file}...")
    cmd = f"bcftools mpileup -f {reference} {bam_file} | bcftools call -mv -Oz -o {output_vcf}"
    subprocess.run(cmd, shell=True, check=True)
    print("Variant calling completed.")

def main():
    # Create output directories
    os.makedirs("results/alignment", exist_ok=True)
    os.makedirs("results/variants", exist_ok=True)
    
    # Reference genome (you'll need to specify the path to your reference genome)
    reference = "data/reference/panthera_tigris.fasta"
    
    # Check if reference genome exists
    if not os.path.exists(reference):
        print(f"Error: Reference genome not found at {reference}")
        print("Please download a reference genome for Panthera tigris and place it at the specified location.")
        return
    
    # Index the reference genome first
    run_bwa_index(reference)
    
    # Process Caspian tiger
    caspian_reads1 = "results/preprocessing/caspian_1_trimmed.fastq.gz"
    caspian_reads2 = "results/preprocessing/caspian_2_trimmed.fastq.gz"
    caspian_bam = "results/alignment/caspian.sorted.bam"
    caspian_vcf = "results/variants/caspian.vcf.gz"
    
    if os.path.exists(caspian_reads1) and os.path.exists(caspian_reads2):
        print("Processing Caspian tiger genome...")
        run_bwa_mem(reference, caspian_reads1, caspian_reads2, caspian_bam)
        run_samtools_index(caspian_bam)
        run_bcftools_call(caspian_bam, caspian_vcf, reference)
    else:
        print(f"Error: Caspian tiger trimmed reads not found")
    
    # Process Amur tiger
    amur_reads1 = "results/preprocessing/amur_1_trimmed.fastq.gz"
    amur_reads2 = "results/preprocessing/amur_2_trimmed.fastq.gz"
    amur_bam = "results/alignment/amur.sorted.bam"
    amur_vcf = "results/variants/amur.vcf.gz"
    
    if os.path.exists(amur_reads1) and os.path.exists(amur_reads2):
        print("Processing Amur tiger genome...")
        run_bwa_mem(reference, amur_reads1, amur_reads2, amur_bam)
        run_samtools_index(amur_bam)
        run_bcftools_call(amur_bam, amur_vcf, reference)
    else:
        print(f"Error: Amur tiger trimmed reads not found")

if __name__ == "__main__":
    main() 