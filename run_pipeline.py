#!/usr/bin/env python3
import os
import subprocess
import sys
import time
from datetime import timedelta

def format_time(seconds):
    """Format seconds into human readable time"""
    return str(timedelta(seconds=int(seconds)))

def run_script(script_name, description):
    """Run a Python script and measure execution time"""
    start_time = time.time()
    print(f"\n{'='*80}")
    print(f"Running {description}...")
    print(f"{'='*80}")
    
    try:
        subprocess.run([sys.executable, script_name], check=True)
        elapsed_time = time.time() - start_time
        print(f"\n{description} completed in {format_time(elapsed_time)}")
        return True
    except subprocess.CalledProcessError as e:
        print(f"Error running {script_name}: {e}")
        return False

def main():
    total_start_time = time.time()
    
    # Create necessary directories
    os.makedirs("data/reference", exist_ok=True)
    os.makedirs("results", exist_ok=True)
    
    # Step 1: Download reference genome
    if not os.path.exists("data/reference/panthera_tigris.fasta"):
        if not run_script("scripts/download_reference.py", "Reference genome download"):
            print("Failed to download reference genome. Exiting.")
            return
    else:
        print("Reference genome already exists. Skipping download.")
    
    # Step 2: Download tiger genomes (if not already downloaded)
    if not (os.path.exists("data/caspian/SRR18572400_1.fastq") and 
            os.path.exists("data/amur/SRR31485304_1.fastq")):
        if not run_script("scripts/download_genomes.py", "Tiger genome download"):
            print("Failed to download tiger genomes. Exiting.")
            return
    else:
        print("Tiger genomes already exist. Skipping download.")
    
    # Step 3: Preprocess the genomes
    if not run_script("scripts/preprocess.py", "Genome preprocessing"):
        print("Failed to preprocess genomes. Exiting.")
        return
    
    # Step 4: Align and call variants
    if not run_script("scripts/align_and_call.py", "Genome alignment and variant calling"):
        print("Failed to align genomes and call variants. Exiting.")
        return
    
    # Step 5: Analyze and visualize results
    if not run_script("scripts/analyze_and_visualize.py", "Comparative analysis and visualization"):
        print("Failed to analyze and visualize results. Exiting.")
        return
    
    total_time = time.time() - total_start_time
    print(f"\n{'='*80}")
    print(f"Pipeline completed successfully in {format_time(total_time)}")
    print(f"{'='*80}")
    print("\nResults are available in the following directories:")
    print("- Preprocessing results: results/preprocessing/")
    print("- Alignment results: results/alignment/")
    print("- Variant calling results: results/variants/")
    print("- Analysis results: results/analysis/")
    print("- Visualizations: results/figures/")

if __name__ == "__main__":
    main() 