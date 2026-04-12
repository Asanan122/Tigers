import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from cyvcf2 import VCF

def load_vcf(vcf_file):
    """Load VCF file and return a pandas DataFrame with variant information"""
    vcf_reader = VCF(vcf_file)
    variants = []
    
    for record in vcf_reader:
        variant = {
            'CHROM': record.CHROM,
            'POS': record.POS,
            'REF': record.REF,
            'ALT': str(record.ALT[0]),
            'QUAL': record.QUAL,
            'FILTER': record.FILTER if record.FILTER else 'PASS',
            'CASPIAN': record.genotypes[0][0],  # First sample (Caspian)
            'AMUR': record.genotypes[1][0]       # Second sample (Amur)
        }
        variants.append(variant)
    
    return pd.DataFrame(variants)

def analyze_variants(df):
    """Analyze variants and return statistics"""
    stats = {
        'total_variants': len(df),
        'shared_variants': len(df[(df['CASPIAN'] != -1) & (df['AMUR'] != -1)]),
        'caspian_unique': len(df[(df['CASPIAN'] != -1) & (df['AMUR'] == -1)]),
        'amur_unique': len(df[(df['CASPIAN'] == -1) & (df['AMUR'] != -1)]),
        'homozygous_differences': len(df[
            ((df['CASPIAN'] == 0) & (df['AMUR'] == 2)) |
            ((df['CASPIAN'] == 2) & (df['AMUR'] == 0))
        ])
    }
    return stats

def plot_variant_distribution(stats):
    """Create pie chart of variant distribution"""
    plt.figure(figsize=(10, 8))
    labels = ['Shared Variants', 'Caspian Unique', 'Amur Unique']
    sizes = [stats['shared_variants'], stats['caspian_unique'], stats['amur_unique']]
    colors = ['#ff9999', '#66b3ff', '#99ff99']
    
    plt.pie(sizes, labels=labels, colors=colors, autopct='%1.1f%%', startangle=90)
    plt.axis('equal')
    plt.title('Distribution of Variants Between Caspian and Amur Tigers')
    plt.savefig('results/analysis/variant_distribution.png')
    plt.close()

def plot_variant_density(df):
    """Create density plot of variants across chromosomes"""
    plt.figure(figsize=(15, 6))
    sns.kdeplot(data=df, x='POS', hue='CHROM', common_norm=False)
    plt.title('Variant Density Across Chromosomes')
    plt.xlabel('Genomic Position')
    plt.ylabel('Density')
    plt.savefig('results/analysis/variant_density.png')
    plt.close()

def main():
    # Create analysis directory if it doesn't exist
    os.makedirs('results/analysis', exist_ok=True)
    
    # Load and analyze variants
    vcf_file = 'results/variants/final.vcf.gz'
    df = load_vcf(vcf_file)
    stats = analyze_variants(df)
    
    # Generate plots
    plot_variant_distribution(stats)
    plot_variant_density(df)
    
    # Save statistics to file
    with open('results/analysis/variant_statistics.txt', 'w') as f:
        f.write("Variant Analysis Statistics\n")
        f.write("=========================\n\n")
        f.write(f"Total Variants: {stats['total_variants']:,}\n")
        f.write(f"Shared Variants: {stats['shared_variants']:,}\n")
        f.write(f"Caspian Unique Variants: {stats['caspian_unique']:,}\n")
        f.write(f"Amur Unique Variants: {stats['amur_unique']:,}\n")
        f.write(f"Homozygous Differences: {stats['homozygous_differences']:,}\n")
    
    print("Analysis completed! Check the results/analysis directory for results.")

if __name__ == "__main__":
    main() 