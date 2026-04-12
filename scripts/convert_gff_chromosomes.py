#!/usr/bin/env python3
import os
import gzip

def load_fasta_headers(fasta_file):
    """Extract chromosome names from FASTA file."""
    headers = set()
    with open(fasta_file, 'r') as f:
        for line in f:
            if line.startswith('>'):
                # Remove '>' and get the first part of the header (chromosome ID)
                chrom = line[1:].split()[0].strip()
                headers.add(chrom)
    return headers

def convert_gff(input_gff, output_gff, ref_chroms):
    """Convert GFF file using reference chromosome names."""
    with open(input_gff, 'r') as infile, open(output_gff, 'w') as outfile:
        for line in infile:
            if line.startswith('#'):
                outfile.write(line)
                continue
            
            fields = line.strip().split('\t')
            if len(fields) < 9:
                outfile.write(line)
                continue
                
            chrom = fields[0]
            # If the chromosome is not in reference, try to find a matching one
            if chrom not in ref_chroms:
                # Try different formats (with/without version)
                chrom_base = chrom.split('.')[0]
                found = False
                for ref_chrom in ref_chroms:
                    if ref_chrom.startswith(chrom_base):
                        fields[0] = ref_chrom
                        found = True
                        break
                if not found:
                    # Skip features with unknown chromosomes
                    continue
            
            outfile.write('\t'.join(fields) + '\n')

def main():
    # Input and output paths
    fasta_file = "data/reference/panthera_tigris.fasta"
    input_gff = "data/reference/panthera_tigris.gff"
    output_gff = "tools/snpEff/data/panthera_tigris/genes.gff"
    
    # Create output directory if it doesn't exist
    os.makedirs(os.path.dirname(output_gff), exist_ok=True)
    
    # Load reference chromosome names
    print("Loading reference chromosome names...")
    ref_chroms = load_fasta_headers(fasta_file)
    print(f"Found {len(ref_chroms)} chromosomes in reference")
    
    # Convert GFF file
    print("Converting GFF file...")
    convert_gff(input_gff, output_gff, ref_chroms)
    print(f"GFF file converted and saved to {output_gff}")

if __name__ == "__main__":
    main() 