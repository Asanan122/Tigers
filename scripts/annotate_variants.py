import os
import subprocess
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import vcfpy
import gzip
import shutil
import pysam
import re

def check_snpeff_installed():
    """Check if SnpEff is installed and available in the PATH"""
    try:
        subprocess.run(['snpeff', '-version'], stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True)
        return True
    except (subprocess.SubprocessError, FileNotFoundError):
        return False

def download_snpeff():
    """Download and install SnpEff"""
    print("SnpEff not found. Downloading and installing...")
    
    # Create a directory for SnpEff
    os.makedirs('tools', exist_ok=True)
    
    # Download SnpEff
    subprocess.run([
        'curl', '-L', 
        'https://sourceforge.net/projects/snpeff/files/snpEff_latest_core.zip/download',
        '-o', 'tools/snpEff.zip'
    ], check=True)
    
    # Unzip SnpEff
    subprocess.run(['unzip', '-o', 'tools/snpEff.zip', '-d', 'tools'], check=True)
    
    # Make SnpEff executable
    subprocess.run(['chmod', '+x', 'tools/snpEff/snpEff.jar'], check=True)
    
    print("SnpEff installed successfully.")
    return 'tools/snpEff/snpEff.jar'

def download_tiger_genome():
    """Download tiger genome for SnpEff database"""
    print("Downloading tiger genome for SnpEff database...")
    
    # Create a directory for the genome
    os.makedirs('data/reference/snpeff_db', exist_ok=True)
    
    # Use the existing reference genome
    genome_file = 'data/reference/panthera_tigris.fasta'
    
    # Create a simple GFF file for annotation
    gff_file = 'data/reference/snpeff_db/panthera_tigris.gff'
    
    # Create a basic GFF file with chromosome information
    with open(gff_file, 'w') as f:
        f.write("##gff-version 3\n")
        f.write("##sequence-region NC_018728.1 1 1000000\n")
        f.write("NC_018728.1\tRefSeq\tchromosome\t1\t1000000\t.\t.\t.\tID=NC_018728.1;Name=NC_018728.1\n")
    
    print("Using existing reference genome and created basic GFF file.")
    return genome_file, gff_file

def create_simplified_gff(input_gff, output_gff):
    """Create simplified GFF with chromosome numbers instead of scaffold IDs."""
    # Read chromosome mapping
    scaffold_to_chr = {}
    with open('data/reference/chromosomes.map', 'r') as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            parts = line.strip().split('\t')
            if len(parts) >= 2:
                scaffold = parts[0]
                chr_num = parts[1]
                scaffold_to_chr[scaffold] = chr_num
    
    with open(input_gff, 'r') as f_in, open(output_gff, 'w') as f_out:
        for line in f_in:
            if line.startswith('#'):
                if line.startswith('##sequence-region'):
                    # Update sequence-region directives
                    parts = line.strip().split()
                    scaffold = parts[1]
                    if scaffold in scaffold_to_chr:
                        chr_num = scaffold_to_chr[scaffold]
                        f_out.write(f'##sequence-region {chr_num} {parts[2]} {parts[3]}\n')
                    else:
                        # If scaffold not in mapping, use the original scaffold ID
                        f_out.write(line)
                else:
                    f_out.write(line)
                continue
            
            parts = line.strip().split('\t')
            if len(parts) >= 9:
                scaffold = parts[0]
                if scaffold in scaffold_to_chr:
                    chr_num = scaffold_to_chr[scaffold]
                    parts[0] = chr_num
                    # Update any ID or Parent attributes that contain scaffold IDs
                    attrs = dict(item.split('=') for item in parts[8].split(';'))
                    if 'ID' in attrs and scaffold in attrs['ID']:
                        attrs['ID'] = attrs['ID'].replace(scaffold, chr_num)
                    if 'Parent' in attrs and scaffold in attrs['Parent']:
                        attrs['Parent'] = attrs['Parent'].replace(scaffold, chr_num)
                    if 'Name' in attrs and scaffold in attrs['Name']:
                        attrs['Name'] = attrs['Name'].replace(scaffold, chr_num)
                    parts[8] = ';'.join(f'{k}={v}' for k, v in attrs.items())
                    f_out.write('\t'.join(parts) + '\n')
                else:
                    # If scaffold not in mapping, write the original line
                    f_out.write(line)

def setup_snpeff_database(genome_file, gff_file):
    """Set up SnpEff database for variant annotation."""
    # Create SnpEff directory structure
    os.makedirs('tools/snpEff/data/panthera_tigris', exist_ok=True)
    
    # Read chromosome mapping
    scaffold_to_chr = {}
    chr_to_scaffold = {}
    with open('data/reference/chromosomes.map', 'r') as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            parts = line.strip().split('\t')
            if len(parts) >= 2:
                scaffold = parts[0]
                chr_num = parts[1]
                scaffold_to_chr[scaffold] = chr_num
                chr_to_scaffold[chr_num] = scaffold
    
    # Create simplified GFF with scaffold IDs
    simplified_gff = 'tools/snpEff/data/panthera_tigris/genes.gff'
    with open(gff_file, 'r') as f_in, open(simplified_gff, 'w') as f_out:
        for line in f_in:
            if line.startswith('#'):
                if line.startswith('##sequence-region'):
                    # Update sequence-region directives
                    parts = line.strip().split()
                    chr_num = parts[1]
                    if chr_num in chr_to_scaffold:
                        scaffold = chr_to_scaffold[chr_num]
                        f_out.write(f'##sequence-region {scaffold} {parts[2]} {parts[3]}\n')
                    continue
                f_out.write(line)
                continue
            
            parts = line.strip().split('\t')
            if len(parts) >= 9:
                chr_num = parts[0]
                if chr_num in chr_to_scaffold:
                    scaffold = chr_to_scaffold[chr_num]
                    parts[0] = scaffold
                    # Update any ID or Parent attributes that contain chromosome numbers
                    attrs = dict(item.split('=') for item in parts[8].split(';'))
                    if 'ID' in attrs and chr_num in attrs['ID']:
                        attrs['ID'] = attrs['ID'].replace(chr_num, scaffold)
                    if 'Parent' in attrs and chr_num in attrs['Parent']:
                        attrs['Parent'] = attrs['Parent'].replace(chr_num, scaffold)
                    if 'Name' in attrs and chr_num in attrs['Name']:
                        attrs['Name'] = attrs['Name'].replace(chr_num, scaffold)
                    parts[8] = ';'.join(f'{k}={v}' for k, v in attrs.items())
                    f_out.write('\t'.join(parts) + '\n')
    
    # Create FASTA file with scaffold IDs
    fasta = pysam.FastaFile(genome_file)
    with open('tools/snpEff/data/panthera_tigris/sequences.fa', 'w') as out:
        for ref in fasta.references:
            out.write(f'>{ref}\n{fasta.fetch(ref)}\n')
    
    # Copy chromosome mapping file to SnpEff directory
    shutil.copy('data/reference/chromosomes.map', 
                'tools/snpEff/data/panthera_tigris/chromosomes.map')
    
    # Update config file
    with open('tools/snpEff/snpEff.config', 'a') as f:
        f.write('\n# Tiger genome, version 1\n')
        f.write('panthera_tigris.genome : Panthera tigris\n')
        f.write('panthera_tigris.chromosomes.map = chromosomes.map\n')
    
    # Build database
    subprocess.run([
        'java', '-jar', 'tools/snpEff/snpEff.jar',
        'build', '-gff3', '-v',
        'panthera_tigris',
        '-c', 'tools/snpEff/snpEff.config'
    ], check=True)

def annotate_variants(input_vcf, output_vcf):
    """Annotate variants using SnpEff."""
    # Create output directory if it doesn't exist
    os.makedirs(os.path.dirname(output_vcf), exist_ok=True)
    
    subprocess.run([
        'java', '-jar', 'tools/snpEff/snpEff.jar',
        'panthera_tigris',
        '-c', 'tools/snpEff/snpEff.config',
        '-chrom', 'chromosomes.map',
        input_vcf
    ], stdout=open(output_vcf, 'w'), check=True)
    
    return output_vcf

def parse_annotated_vcf(vcf_file):
    """Parse the annotated VCF file and extract variant information"""
    print("Parsing annotated VCF file...")
    
    # Create a temporary file without the header
    temp_file = 'results/variants/annotated/temp.vcf'
    
    with gzip.open(vcf_file, 'rt') as f_in, open(temp_file, 'w') as f_out:
        for line in f_in:
            if not line.startswith('#'):
                f_out.write(line)
    
    # Read the VCF file
    reader = vcfpy.Reader.from_path(temp_file)
    records = []
    
    for record in reader:
        # Get genotype information
        caspian_call = record.call_for_sample['caspian']
        amur_call = record.call_for_sample['amur']
        
        # Extract annotation information
        ann = record.INFO.get('ANN', [''])[0].split('|') if 'ANN' in record.INFO else [''] * 16
        
        # Extract gene information
        gene = ann[3] if len(ann) > 3 else 'Unknown'
        gene_id = ann[4] if len(ann) > 4 else 'Unknown'
        feature_type = ann[5] if len(ann) > 5 else 'Unknown'
        impact = ann[2] if len(ann) > 2 else 'Unknown'
        effect = ann[1] if len(ann) > 1 else 'Unknown'
        hgvs_c = ann[10] if len(ann) > 10 else 'Unknown'
        hgvs_p = ann[11] if len(ann) > 11 else 'Unknown'
        
        records.append({
            'CHROM': record.CHROM,
            'POS': record.POS,
            'REF': record.REF,
            'ALT': str(record.ALT[0].value),
            'QUAL': record.QUAL if record.QUAL is not None else 0,
            'FILTER': ','.join(record.FILTER) if record.FILTER else 'PASS',
            'gene': gene,
            'gene_id': gene_id,
            'feature_type': feature_type,
            'impact': impact,
            'effect': effect,
            'hgvs_c': hgvs_c,
            'hgvs_p': hgvs_p,
            'caspian_gt': caspian_call.gt_type,
            'amur_gt': amur_call.gt_type,
            'caspian_gq': caspian_call.data.get('GQ', 0),
            'amur_gq': amur_call.data.get('GQ', 0)
        })
    
    # Remove temporary file
    os.remove(temp_file)
    
    return pd.DataFrame(records)

def analyze_annotated_variants(df):
    """Analyze the annotated variants"""
    print("Analyzing annotated variants...")
    
    # Create results directory
    Path('results/analysis/annotated').mkdir(parents=True, exist_ok=True)
    
    # Basic variant statistics
    stats = {
        'total_variants': len(df),
        'variants_in_genes': len(df[df['gene'] != 'Unknown']),
        'high_impact_variants': len(df[df['impact'] == 'HIGH']),
        'moderate_impact_variants': len(df[df['impact'] == 'MODERATE']),
        'low_impact_variants': len(df[df['impact'] == 'LOW']),
        'modifier_impact_variants': len(df[df['impact'] == 'MODIFIER']),
        'shared_variants': len(df[df['caspian_gt'] == df['amur_gt']]),
        'caspian_unique': len(df[df['caspian_gt'] != df['amur_gt']]),
        'amur_unique': len(df[df['amur_gt'] != df['caspian_gt']])
    }
    
    # Generate detailed CSV files
    generate_annotated_reports(df)
    
    # Generate visualizations
    generate_annotated_visualizations(df)
    
    return stats

def generate_annotated_reports(df):
    """Generate detailed CSV reports for annotated variants"""
    # 1. All annotated variants
    df.to_csv('results/analysis/annotated/all_annotated_variants.csv', index=False)
    
    # 2. High-impact variants
    high_impact = df[df['impact'] == 'HIGH']
    high_impact.to_csv('results/analysis/annotated/high_impact_variants.csv', index=False)
    
    # 3. Variants in genes
    gene_variants = df[df['gene'] != 'Unknown']
    gene_variants.to_csv('results/analysis/annotated/gene_variants.csv', index=False)
    
    # 4. Unique variants by subspecies
    caspian_unique = df[df['caspian_gt'] != df['amur_gt']]
    caspian_unique.to_csv('results/analysis/annotated/caspian_unique_variants.csv', index=False)
    
    amur_unique = df[df['amur_gt'] != df['caspian_gt']]
    amur_unique.to_csv('results/analysis/annotated/amur_unique_variants.csv', index=False)
    
    # 5. Gene-wise summary
    gene_summary = df[df['gene'] != 'Unknown'].groupby('gene').agg({
        'impact': 'value_counts',
        'effect': 'value_counts',
        'POS': 'count'
    }).reset_index()
    gene_summary.to_csv('results/analysis/annotated/gene_summary.csv', index=False)

def generate_annotated_visualizations(df):
    """Generate visualizations for annotated variants"""
    # 1. Impact distribution
    plt.figure(figsize=(10, 6))
    df['impact'].value_counts().plot(kind='bar')
    plt.title('Distribution of Variant Impact')
    plt.xlabel('Impact')
    plt.ylabel('Count')
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.savefig('results/analysis/annotated/impact_distribution.png')
    plt.close()
    
    # 2. Effect distribution
    plt.figure(figsize=(12, 6))
    df['effect'].value_counts().head(10).plot(kind='bar')
    plt.title('Top 10 Most Common Variant Effects')
    plt.xlabel('Effect')
    plt.ylabel('Count')
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.savefig('results/analysis/annotated/effect_distribution.png')
    plt.close()
    
    # 3. Feature type distribution
    plt.figure(figsize=(10, 6))
    df['feature_type'].value_counts().plot(kind='bar')
    plt.title('Distribution of Feature Types')
    plt.xlabel('Feature Type')
    plt.ylabel('Count')
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.savefig('results/analysis/annotated/feature_type_distribution.png')
    plt.close()
    
    # 4. Gene mutation frequency (top 20 genes)
    plt.figure(figsize=(15, 6))
    gene_counts = df[df['gene'] != 'Unknown']['gene'].value_counts().head(20)
    gene_counts.plot(kind='bar')
    plt.title('Top 20 Genes with Most Variants')
    plt.xlabel('Gene')
    plt.ylabel('Number of Variants')
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.savefig('results/analysis/annotated/gene_frequency.png')
    plt.close()
    
    # 5. Impact by feature type
    plt.figure(figsize=(12, 6))
    impact_by_feature = pd.crosstab(df['feature_type'], df['impact'])
    impact_by_feature.plot(kind='bar', stacked=True)
    plt.title('Impact Distribution by Feature Type')
    plt.xlabel('Feature Type')
    plt.ylabel('Count')
    plt.xticks(rotation=45)
    plt.legend(title='Impact')
    plt.tight_layout()
    plt.savefig('results/analysis/annotated/impact_by_feature.png')
    plt.close()

def main():
    # Check if SnpEff is installed
    snpeff_jar = 'snpeff' if check_snpeff_installed() else download_snpeff()
    
    # Download tiger genome if needed
    genome_file = 'data/reference/panthera_tigris.fasta'
    gff_file = 'data/reference/snpeff_db/panthera_tigris.gff'
    
    if not os.path.exists(gff_file):
        genome_file, gff_file = download_tiger_genome()
    
    # Set up SnpEff database
    setup_snpeff_database(genome_file, gff_file)
    
    # Annotate variants
    vcf_file = 'results/variants/final.vcf.gz'
    annotated_vcf = 'results/variants/annotated/final.annotated.vcf'
    annotated_vcf = annotate_variants(vcf_file, annotated_vcf)
    
    # Parse annotated VCF
    df = parse_annotated_vcf(annotated_vcf)
    
    # Analyze annotated variants
    stats = analyze_annotated_variants(df)
    
    # Save statistics
    with open('results/analysis/annotated/annotated_statistics.txt', 'w') as f:
        f.write("Annotated Variant Analysis Statistics\n")
        f.write("===================================\n\n")
        for key, value in stats.items():
            f.write(f"{key}: {value:,}\n")
    
    print("Annotation and analysis complete! Results are available in the results/analysis/annotated directory.")
    print("\nGenerated files:")
    print("1. annotated_statistics.txt - Overall statistics")
    print("2. all_annotated_variants.csv - Complete annotated variant information")
    print("3. high_impact_variants.csv - Variants with high impact")
    print("4. gene_variants.csv - Variants in genes")
    print("5. caspian_unique_variants.csv - Variants unique to Caspian tiger")
    print("6. amur_unique_variants.csv - Variants unique to Amur tiger")
    print("7. gene_summary.csv - Gene-wise variant summary")
    print("\nGenerated visualizations:")
    print("1. impact_distribution.png - Distribution of variant impact")
    print("2. effect_distribution.png - Distribution of variant effects")
    print("3. feature_type_distribution.png - Distribution of feature types")
    print("4. gene_frequency.png - Top 20 genes with most variants")
    print("5. impact_by_feature.png - Impact distribution by feature type")

if __name__ == "__main__":
    main() 