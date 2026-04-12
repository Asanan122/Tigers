import os
import subprocess
import sys
import time
from pathlib import Path
from datetime import timedelta

def format_time(seconds):
    """Format seconds into human readable time"""
    return str(timedelta(seconds=int(seconds)))

def download_sra(sra_id, output_dir):
    """
    Download and convert SRA to FASTQ using sratoolkit
    """
    start_time = time.time()
    
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    # Use fasterq-dump with progress information
    cmd = f"fasterq-dump --progress --outdir {output_dir} {sra_id}"
    print(f"\nDownloading and converting {sra_id}...")
    print("This may take a while. Progress will be shown below:")
    subprocess.run(cmd, shell=True, check=True)
    
    # Compress the FASTQ file
    fastq_file = os.path.join(output_dir, f"{sra_id}.fastq")
    if os.path.exists(fastq_file):
        print(f"\nCompressing {fastq_file}...")
        subprocess.run(f"gzip {fastq_file}", shell=True, check=True)
        elapsed_time = time.time() - start_time
        print(f"\nCompleted in {format_time(elapsed_time)}")
        return os.path.join(output_dir, f"{sra_id}.fastq.gz")
    else:
        raise FileNotFoundError(f"Failed to download {sra_id}")

def main():
    total_start_time = time.time()
    
    # Create directories if they don't exist
    os.makedirs('data/caspian', exist_ok=True)
    os.makedirs('data/amur', exist_ok=True)
    
    # Download Caspian tiger genome
    caspian_sra = "SRR18572400"
    try:
        print("\n=== Downloading Caspian Tiger Genome ===")
        caspian_output = download_sra(caspian_sra, "data/caspian")
        print(f"Successfully downloaded Caspian tiger genome to {caspian_output}")
    except Exception as e:
        print(f"Error downloading Caspian tiger genome: {e}")
        sys.exit(1)
    
    # Download Amur tiger genome
    amur_sra = "SRR31485304"
    try:
        print("\n=== Downloading Amur Tiger Genome ===")
        amur_output = download_sra(amur_sra, "data/amur")
        print(f"Successfully downloaded Amur tiger genome to {amur_output}")
    except Exception as e:
        print(f"Error downloading Amur tiger genome: {e}")
        sys.exit(1)
    
    total_time = time.time() - total_start_time
    print(f"\nTotal download time: {format_time(total_time)}")

if __name__ == "__main__":
    main() 