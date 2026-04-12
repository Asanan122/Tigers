import os
import pandas as pd
import numpy as np
import vcfpy
import gzip
from pathlib import Path
import matplotlib.pyplot as plt
import seaborn as sns

def parse_vcf_to_dataframe(vcf_file):
    """Parse VCF file and convert to pandas DataFrame with detailed variant information."""
    print(f"Parsing VCF file: {vcf_file}")
    
    records = []
    with gzip.open(vcf_file, 'rt') as f:
        for line in f:
            if line.startswith('#'):
                continue
                
            parts = line.strip().split('\t')
            chrom = parts[0]
            pos = int(parts[1])
            ref = parts[3]
            alt = parts[4]
            qual = float(parts[5]) if parts[5] != '.' else np.nan
            filter_status = parts[6]
            
            # Parse INFO field
            info = {}
            for item in parts[7].split(';'):
                if '=' in item:
                    key, value = item.split('=', 1)
                    info[key] = value
            
            # Parse genotype information
            format_fields = parts[8].split(':')
            caspian_gt = parts[9].split(':')
            amur_gt = parts[10].split(':')
            
            # Create genotype dictionary
            caspian_genotype = dict(zip(format_fields, caspian_gt))
            amur_genotype = dict(zip(format_fields, amur_gt))
            
            # Extract annotation information if available
            ann = info.get('ANN', '').split('|') if 'ANN' in info else [''] * 16
            gene = ann[3] if len(ann) > 3 else 'Unknown'
            gene_id = ann[4] if len(ann) > 4 else 'Unknown'
            feature_type = ann[5] if len(ann) > 5 else 'Unknown'
            impact = ann[2] if len(ann) > 2 else 'Unknown'
            effect = ann[1] if len(ann) > 1 else 'Unknown'
            
            record = {
                'chromosome': chrom,
                'position': pos,
                'reference': ref,
                'alternative': alt,
                'quality': qual,
                'filter': filter_status,
                'gene': gene,
                'gene_id': gene_id,
                'feature_type': feature_type,
                'impact': impact,
                'effect': effect,
                'caspian_genotype': caspian_genotype.get('GT', './.'),
                'amur_genotype': amur_genotype.get('GT', './.'),
                'caspian_depth': int(caspian_genotype.get('DP', 0)),
                'amur_depth': int(amur_genotype.get('DP', 0)),
                'caspian_quality': float(caspian_genotype.get('GQ', 0)),
                'amur_quality': float(amur_genotype.get('GQ', 0))
            }
            records.append(record)
    
    return pd.DataFrame(records)

def compare_genotypes(df):
    """Compare genotypes between Caspian and Amur tigers and add comparison columns."""
    print("Comparing genotypes between Caspian and Amur tigers...")
    
    # Add comparison columns
    df['genotype_match'] = df['caspian_genotype'] == df['amur_genotype']
    df['variant_type'] = 'shared'
    df.loc[df['caspian_genotype'] == './.', 'variant_type'] = 'amur_unique'
    df.loc[df['amur_genotype'] == './.', 'variant_type'] = 'caspian_unique'
    
    # Add zygosity information
    df['caspian_zygosity'] = df['caspian_genotype'].apply(lambda x: 'homozygous' if x in ['0/0', '1/1'] else 'heterozygous' if x in ['0/1', '1/0'] else 'missing')
    df['amur_zygosity'] = df['amur_genotype'].apply(lambda x: 'homozygous' if x in ['0/0', '1/1'] else 'heterozygous' if x in ['0/1', '1/0'] else 'missing')
    
    return df

def generate_comparison_report(df, output_file):
    """Generate a detailed CSV report of genome comparison."""
    print(f"Generating comparison report: {output_file}")
    
    # Select and order columns for the report
    report_columns = [
        'chromosome', 'position', 'reference', 'alternative',
        'gene', 'gene_id', 'feature_type', 'impact', 'effect',
        'variant_type', 'genotype_match',
        'caspian_genotype', 'caspian_zygosity', 'caspian_depth', 'caspian_quality',
        'amur_genotype', 'amur_zygosity', 'amur_depth', 'amur_quality'
    ]
    
    # Generate the report
    df[report_columns].to_csv(output_file, index=False)
    
    # Print summary statistics
    print("\nComparison Summary:")
    print(f"Total variants: {len(df)}")
    print(f"Shared variants: {len(df[df['variant_type'] == 'shared'])}")
    print(f"Caspian-specific variants: {len(df[df['variant_type'] == 'caspian_unique'])}")
    print(f"Amur-specific variants: {len(df[df['variant_type'] == 'amur_unique'])}")
    
    # Print impact distribution
    print("\nImpact Distribution:")
    print(df['impact'].value_counts())
    
    # Print effect distribution
    print("\nEffect Distribution:")
    print(df['effect'].value_counts())

def generate_comparison_plots(df, output_dir):
    """Generate visualization plots for genome comparison."""
    print(f"Generating comparison plots in: {output_dir}")
    os.makedirs(output_dir, exist_ok=True)
    
    # Set style
    plt.style.use('seaborn')
    
    # 1. Variant Type Distribution
    plt.figure(figsize=(10, 6))
    df['variant_type'].value_counts().plot(kind='bar')
    plt.title('Distribution of Variant Types')
    plt.xlabel('Variant Type')
    plt.ylabel('Count')
    plt.savefig(os.path.join(output_dir, 'variant_type_distribution.png'))
    plt.close()
    
    # 2. Impact Distribution
    plt.figure(figsize=(12, 6))
    df['impact'].value_counts().plot(kind='bar')
    plt.title('Distribution of Variant Impacts')
    plt.xlabel('Impact')
    plt.ylabel('Count')
    plt.xticks(rotation=45)
    plt.savefig(os.path.join(output_dir, 'impact_distribution.png'))
    plt.close()
    
    # 3. Depth Comparison
    plt.figure(figsize=(10, 6))
    plt.scatter(df['caspian_depth'], df['amur_depth'], alpha=0.5)
    plt.plot([0, max(df['caspian_depth'])], [0, max(df['amur_depth'])], 'r--')
    plt.title('Read Depth Comparison')
    plt.xlabel('Caspian Depth')
    plt.ylabel('Amur Depth')
    plt.savefig(os.path.join(output_dir, 'depth_comparison.png'))
    plt.close()

def main():
    # Create output directories
    os.makedirs('results/comparison', exist_ok=True)
    os.makedirs('results/comparison/plots', exist_ok=True)
    
    # Input VCF file
    vcf_file = 'results/variants/annotated/annotated.vcf.gz'
    
    # Parse VCF to DataFrame
    df = parse_vcf_to_dataframe(vcf_file)
    
    # Compare genotypes
    df = compare_genotypes(df)
    
    # Generate comparison report
    output_csv = 'results/comparison/genome_comparison.csv'
    generate_comparison_report(df, output_csv)
    
    # Generate comparison plots
    generate_comparison_plots(df, 'results/comparison/plots')
    
    print("\nGenome comparison completed successfully!")
    print(f"Detailed comparison report saved to: {output_csv}")
    print("Visualization plots saved to: results/comparison/plots/")

if __name__ == "__main__":
    main() 