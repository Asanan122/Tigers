import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import vcfpy
import os
from Bio import SeqIO
from collections import defaultdict

def read_vcf(vcf_file):
    """Read VCF file and return a pandas DataFrame with detailed variant information"""
    reader = vcfpy.Reader.from_path(vcf_file)
    records = []
    
    for record in reader:
        # Get genotype information
        caspian_call = record.call_for_sample['caspian']
        amur_call = record.call_for_sample['amur']
        
        # Calculate variant effect
        variant_type = get_variant_type(record.REF, str(record.ALT[0].value))
        
        # Extract variant impact from INFO field if available
        impact = record.INFO.get('ANN', [''])[0].split('|')[2] if 'ANN' in record.INFO else 'Unknown'
        gene = record.INFO.get('ANN', [''])[0].split('|')[3] if 'ANN' in record.INFO else 'Unknown'
        
        records.append({
            'CHROM': record.CHROM,
            'POS': record.POS,
            'REF': record.REF,
            'ALT': str(record.ALT[0].value),
            'QUAL': record.QUAL if record.QUAL is not None else 0,
            'FILTER': ','.join(record.FILTER) if record.FILTER else 'PASS',
            'variant_type': variant_type,
            'impact': impact,
            'gene': gene,
            'caspian_gt': caspian_call.gt_type,
            'amur_gt': amur_call.gt_type,
            'caspian_gq': caspian_call.data.get('GQ', 0),
            'amur_gq': amur_call.data.get('GQ', 0)
        })
    
    return pd.DataFrame(records)

def get_variant_type(ref, alt):
    """Determine the type of variant"""
    if len(ref) == len(alt) == 1:
        return 'SNP'
    elif len(ref) > len(alt):
        return 'Deletion'
    elif len(ref) < len(alt):
        return 'Insertion'
    else:
        return 'Complex'

def analyze_variants(df):
    """Perform detailed variant analysis"""
    # Create results directory
    Path('results/analysis').mkdir(parents=True, exist_ok=True)
    
    # Basic variant statistics
    stats = {
        'total_variants': len(df),
        'snps': len(df[df['variant_type'] == 'SNP']),
        'insertions': len(df[df['variant_type'] == 'Insertion']),
        'deletions': len(df[df['variant_type'] == 'Deletion']),
        'shared_variants': len(df[df['caspian_gt'] == df['amur_gt']]),
        'caspian_unique': len(df[df['caspian_gt'] != df['amur_gt']]),
        'amur_unique': len(df[df['amur_gt'] != df['caspian_gt']]),
        'high_impact_variants': len(df[df['impact'] == 'HIGH'])
    }
    
    # Generate detailed CSV files
    generate_detailed_reports(df)
    
    # Generate visualizations
    generate_visualizations(df)
    
    return stats

def generate_detailed_reports(df):
    """Generate detailed CSV reports for different aspects of the analysis"""
    # 1. All variants with details
    df.to_csv('results/analysis/all_variants.csv', index=False)
    
    # 2. High-impact variants
    high_impact = df[df['impact'] == 'HIGH']
    high_impact.to_csv('results/analysis/high_impact_variants.csv', index=False)
    
    # 3. Unique variants by subspecies
    caspian_unique = df[df['caspian_gt'] != df['amur_gt']]
    caspian_unique.to_csv('results/analysis/caspian_unique_variants.csv', index=False)
    
    amur_unique = df[df['amur_gt'] != df['caspian_gt']]
    amur_unique.to_csv('results/analysis/amur_unique_variants.csv', index=False)
    
    # 4. Gene-wise summary
    gene_summary = df.groupby('gene').agg({
        'variant_type': 'value_counts',
        'impact': 'value_counts',
        'POS': 'count'
    }).reset_index()
    gene_summary.to_csv('results/analysis/gene_summary.csv', index=False)

def generate_visualizations(df):
    """Generate comprehensive visualizations of the variant data"""
    # 1. Variant type distribution
    plt.figure(figsize=(10, 6))
    df['variant_type'].value_counts().plot(kind='bar')
    plt.title('Distribution of Variant Types')
    plt.xlabel('Variant Type')
    plt.ylabel('Count')
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.savefig('results/analysis/variant_types.png')
    plt.close()
    
    # 2. Variant distribution by chromosome
    plt.figure(figsize=(12, 6))
    df['CHROM'].value_counts().plot(kind='bar')
    plt.title('Variant Distribution by Chromosome')
    plt.xlabel('Chromosome')
    plt.ylabel('Number of Variants')
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.savefig('results/analysis/chromosome_distribution.png')
    plt.close()
    
    # 3. Quality score distribution
    plt.figure(figsize=(10, 6))
    sns.histplot(data=df, x='QUAL', bins=50)
    plt.title('Distribution of Variant Quality Scores')
    plt.xlabel('Quality Score')
    plt.ylabel('Count')
    plt.tight_layout()
    plt.savefig('results/analysis/quality_distribution.png')
    plt.close()
    
    # 4. Impact distribution
    plt.figure(figsize=(10, 6))
    df['impact'].value_counts().plot(kind='bar')
    plt.title('Distribution of Variant Impact')
    plt.xlabel('Impact')
    plt.ylabel('Count')
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.savefig('results/analysis/impact_distribution.png')
    plt.close()
    
    # 5. Gene mutation frequency (top 20 genes)
    plt.figure(figsize=(15, 6))
    gene_counts = df['gene'].value_counts().head(20)
    gene_counts.plot(kind='bar')
    plt.title('Top 20 Genes with Most Variants')
    plt.xlabel('Gene')
    plt.ylabel('Number of Variants')
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.savefig('results/analysis/gene_frequency.png')
    plt.close()

def main():
    # Read the final VCF file
    vcf_file = 'results/variants/final.vcf.gz'
    if not os.path.exists(vcf_file):
        print(f"Error: {vcf_file} not found")
        return
    
    print("Reading VCF file...")
    df = read_vcf(vcf_file)
    
    print("Analyzing variants...")
    stats = analyze_variants(df)
    
    # Save statistics
    print("Saving results...")
    with open('results/analysis/detailed_statistics.txt', 'w') as f:
        f.write("Detailed Variant Analysis Statistics\n")
        f.write("==================================\n\n")
        for key, value in stats.items():
            f.write(f"{key}: {value:,}\n")
    
    print("Analysis complete! Results are available in the results/analysis directory.")
    print("\nGenerated files:")
    print("1. detailed_statistics.txt - Overall statistics")
    print("2. all_variants.csv - Complete variant information")
    print("3. high_impact_variants.csv - Variants with high impact")
    print("4. caspian_unique_variants.csv - Variants unique to Caspian tiger")
    print("5. amur_unique_variants.csv - Variants unique to Amur tiger")
    print("6. gene_summary.csv - Gene-wise variant summary")
    print("\nGenerated visualizations:")
    print("1. variant_types.png - Distribution of variant types")
    print("2. chromosome_distribution.png - Variants per chromosome")
    print("3. quality_distribution.png - Quality score distribution")
    print("4. impact_distribution.png - Impact level distribution")
    print("5. gene_frequency.png - Top 20 genes with most variants")

if __name__ == "__main__":
    main() 