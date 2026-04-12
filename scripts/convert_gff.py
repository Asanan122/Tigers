#!/usr/bin/env python3
import os
from collections import defaultdict

def load_chromosome_map(map_file):
    """Load chromosome mapping from the map file."""
    chrom_map = {}
    with open(map_file, 'r') as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            ref_chrom, gff_chrom = line.strip().split()
            chrom_map[ref_chrom] = gff_chrom
    return chrom_map

def convert_gff(input_gff, output_gff, chrom_map):
    """Convert GFF file using chromosome mapping."""
    with open(input_gff, 'r') as infile, open(output_gff, 'w') as outfile:
        for line in infile:
            if line.startswith('#'):
                outfile.write(line)
                continue
            
            fields = line.strip().split('\t')
            if len(fields) < 9:
                outfile.write(line)
                continue
                
            ref_chrom = fields[0]
            if ref_chrom in chrom_map:
                fields[0] = chrom_map[ref_chrom]
            else:
                fields[0] = 'UNKNOWN'
            
            outfile.write('\t'.join(fields) + '\n')

def main():
    # Input and output paths
    input_gff = "data/reference/panthera_tigris.gff"
    output_gff = "tools/snpEff/data/panthera_tigris/genes.gff"
    chrom_map_file = "data/reference/chromosomes.map"
    
    # Create output directory if it doesn't exist
    os.makedirs(os.path.dirname(output_gff), exist_ok=True)
    
    # Load chromosome mapping
    chrom_map = load_chromosome_map(chrom_map_file)
    
    # Convert GFF file
    convert_gff(input_gff, output_gff, chrom_map)
    print(f"GFF file converted and saved to {output_gff}")

if __name__ == "__main__":
    main() 