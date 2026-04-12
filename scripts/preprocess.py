import os
import subprocess
from pathlib import Path

def run_fastqc(input_file, output_dir):
    """Run FastQC on input file"""
    cmd = f"fastqc {input_file} -o {output_dir}"
    subprocess.run(cmd, shell=True, check=True)

def run_trimmomatic(input_file1, input_file2, output_prefix):
    """Run Trimmomatic for adapter trimming and quality filtering on paired-end reads"""
    cmd = f"trimmomatic PE {input_file1} {input_file2} " \
          f"{output_prefix}_1_trimmed.fastq.gz {output_prefix}_1_unpaired.fastq.gz " \
          f"{output_prefix}_2_trimmed.fastq.gz {output_prefix}_2_unpaired.fastq.gz " \
          f"LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36"
    subprocess.run(cmd, shell=True, check=True)

def main():
    # Create output directories
    os.makedirs("results/preprocessing", exist_ok=True)
    
    # Process Caspian tiger genome
    caspian_input1 = "data/caspian/SRR18572400_1.fastq"
    caspian_input2 = "data/caspian/SRR18572400_2.fastq"
    
    if os.path.exists(caspian_input1) and os.path.exists(caspian_input2):
        print("Processing Caspian tiger genome...")
        run_fastqc(caspian_input1, "results/preprocessing")
        run_fastqc(caspian_input2, "results/preprocessing")
        run_trimmomatic(caspian_input1, caspian_input2, "results/preprocessing/caspian")
    else:
        print(f"Error: Caspian tiger input files not found")
    
    # Process Amur tiger genome
    amur_input1 = "data/amur/SRR31485304_1.fastq"
    amur_input2 = "data/amur/SRR31485304_2.fastq"
    
    if os.path.exists(amur_input1) and os.path.exists(amur_input2):
        print("Processing Amur tiger genome...")
        run_fastqc(amur_input1, "results/preprocessing")
        run_fastqc(amur_input2, "results/preprocessing")
        run_trimmomatic(amur_input1, amur_input2, "results/preprocessing/amur")
    else:
        print(f"Error: Amur tiger input files not found")

if __name__ == "__main__":
    main() 