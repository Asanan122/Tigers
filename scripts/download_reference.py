import os
import subprocess
import sys
import time
from pathlib import Path
from datetime import timedelta

def format_time(seconds):
    """Format seconds into human readable time"""
    return str(timedelta(seconds=int(seconds)))

def download_reference(output_dir):
    """
    Download a reference genome for Panthera tigris
    """
    start_time = time.time()
    
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    # Using the Amur Tiger reference genome from Ensembl
    output_file = os.path.join(output_dir, "panthera_tigris.fasta")
    compressed_file = f"{output_file}.gz"
    
    print(f"\nDownloading reference genome for Panthera tigris...")
    print("This may take a while. Progress will be shown below:")
    
    # Use curl to download the reference genome from Ensembl
    cmd = f"curl -L -o {compressed_file} http://ftp.ensembl.org/pub/release-110/fasta/panthera_tigris_altaica/dna/Panthera_tigris_altaica.PanTig1.0.dna.toplevel.fa.gz"
    subprocess.run(cmd, shell=True, check=True)
    
    # Decompress the file
    print(f"\nDecompressing reference genome...")
    cmd = f"gunzip -c {compressed_file} > {output_file}"
    subprocess.run(cmd, shell=True, check=True)
    
    # Remove the compressed file
    os.remove(compressed_file)
    
    elapsed_time = time.time() - start_time
    print(f"\nCompleted in {format_time(elapsed_time)}")
    
    return output_file

def main():
    # Create reference directory
    ref_dir = "data/reference"
    os.makedirs(ref_dir, exist_ok=True)
    
    try:
        ref_file = download_reference(ref_dir)
        print(f"Successfully downloaded reference genome to {ref_file}")
    except Exception as e:
        print(f"Error downloading reference genome: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main() 