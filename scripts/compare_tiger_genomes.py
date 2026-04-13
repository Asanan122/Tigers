import pandas as pd
import os
from collections import defaultdict
import matplotlib.pyplot as plt
import seaborn as sns
from tqdm import tqdm
import cyvcf2
import logging
import numpy as np
from typing import Dict, List, Tuple
import gzip
import glob
import sys

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler('variant_analysis.log'),
        logging.StreamHandler()
    ]
)

def parse_gff_file(gff_file: str) -> pd.DataFrame:
    """Parse GFF file directly and extract gene information."""
    genes = []
    logging.info("Parsing GFF file...")
    try:
        with open(gff_file, 'r') as f:
            for line in tqdm(f, desc="Reading GFF"):
                if line.startswith('#'):
                    continue
                fields = line.strip().split('\t')
                if len(fields) < 9 or fields[2] != 'gene':
                    continue
                    
                # Parse attributes
                attrs = dict(item.split('=') for item in fields[8].split(';') if '=' in item)
                
                genes.append({
                    'gene_id': attrs.get('ID', ''),
                    'chromosome': fields[0],
                    'start': int(fields[3]),
                    'end': int(fields[4]),
                    'strand': fields[6],
                    'name': attrs.get('Name', ''),
                    'description': attrs.get('description', '')
                })
        logging.info(f"Successfully parsed {len(genes)} genes from GFF file")
        return pd.DataFrame(genes)
    except Exception as e:
        logging.error(f"Error parsing GFF file: {str(e)}")
        raise

def get_variants(vcf_file: str, chromosomes: List[str]) -> Dict[str, List[dict]]:
    """Extract variants from VCF file using cyvcf2 for efficient parsing."""
    variants = defaultdict(list)
    logging.info(f"Processing VCF file: {vcf_file}")
    
    try:
        vcf = cyvcf2.VCF(vcf_file)
        # Get all chromosomes from VCF
        vcf_chroms = set(vcf.seqnames)
        logging.info(f"Found {len(vcf_chroms)} chromosomes in VCF file")
        
        # Process one chromosome at a time
        for chrom in tqdm(vcf_chroms, desc="Processing chromosomes"):
            chrom_variants = []
            variant_count = 0
            
            try:
                # Get variants for this chromosome
                for variant in vcf(chrom):
                    # Get position and alleles
                    pos = variant.POS
                    ref = variant.REF
                    alt = variant.ALT[0] if variant.ALT else None
                    
                    if not alt:
                        continue
                        
                    # Get quality
                    quality = variant.QUAL if variant.QUAL else 0
                    
                    # Get depth from INFO
                    depth = variant.INFO.get('DP', 0)
                    
                    # Determine variant type
                    variant_type = 'SNP' if len(ref) == len(alt) == 1 else 'INDEL'
                    
                    chrom_variants.append({
                        'position': pos,
                        'ref': ref,
                        'alt': alt,
                        'quality': quality,
                        'depth': depth,
                        'type': variant_type
                    })
                    
                    variant_count += 1
                    
                    # Process in chunks to manage memory
                    if variant_count >= 10000:  # Smaller chunk size
                        variants[chrom].extend(chrom_variants)
                        chrom_variants = []
                        variant_count = 0
                        # Force garbage collection
                        import gc
                        gc.collect()
                
                # Add remaining variants
                if chrom_variants:
                    variants[chrom].extend(chrom_variants)
                
                logging.info(f"Processed {len(variants[chrom])} variants for chromosome {chrom}")
                
            except Exception as e:
                logging.warning(f"Warning processing chromosome {chrom}: {str(e)}")
                continue
            
        vcf.close()
        total_variants = sum(len(v) for v in variants.values())
        logging.info(f"Successfully processed {total_variants} variants across {len(variants)} chromosomes")
        return variants
    except Exception as e:
        logging.error(f"Error processing VCF file: {str(e)}")
        raise

def load_chromosome_mapping(mapping_file: str) -> Dict[str, str]:
    """Load chromosome mapping from file."""
    mapping = {}
    try:
        with open(mapping_file, 'r') as f:
            for line in f:
                chrom = line.strip()
                if chrom:
                    # Map both ways
                    mapping[chrom] = chrom
        logging.info(f"Loaded {len(mapping)} chromosome mappings")
        return mapping
    except Exception as e:
        logging.error(f"Error loading chromosome mapping: {str(e)}")
        raise

def analyze_gene(gene: pd.Series, caspian_variants: Dict[str, List[dict]], 
                amur_variants: Dict[str, List[dict]], chrom_mapping: Dict[str, str]) -> dict:
    """Analyze variants within a gene region."""
    try:
        chrom = gene['chromosome']
        start = gene['start']
        end = gene['end']
        
        # Initialize counters
        caspian_snps = 0
        caspian_indels = 0
        amur_snps = 0
        amur_indels = 0
        caspian_count = 0
        amur_count = 0
        
        # Get mapped chromosome name
        mapped_chrom = chrom_mapping.get(chrom, chrom)
        
        # Process variants directly without creating intermediate lists
        for vcf_chrom, variants in caspian_variants.items():
            if vcf_chrom == mapped_chrom:
                for v in variants:
                    if start <= v['position'] <= end:
                        caspian_count += 1
                        if v['type'] == 'SNP':
                            caspian_snps += 1
                        else:
                            caspian_indels += 1
        
        for vcf_chrom, variants in amur_variants.items():
            if vcf_chrom == mapped_chrom:
                for v in variants:
                    if start <= v['position'] <= end:
                        amur_count += 1
                        if v['type'] == 'SNP':
                            amur_snps += 1
                        else:
                            amur_indels += 1
        
        # Calculate variation metrics
        total_variants = caspian_count + amur_count
        gene_length = end - start + 1
        variation_per_kb = (total_variants / gene_length) * 1000 if gene_length > 0 else 0
        
        return {
            'gene_id': gene['gene_id'],
            'name': gene['name'],
            'description': gene['description'],
            'chromosome': chrom,
            'start': start,
            'end': end,
            'length': gene_length,
            'caspian_variants': caspian_count,
            'amur_variants': amur_count,
            'caspian_snps': caspian_snps,
            'caspian_indels': caspian_indels,
            'amur_snps': amur_snps,
            'amur_indels': amur_indels,
            'total_variants': total_variants,
            'variation_per_kb': variation_per_kb
        }
        
    except Exception as e:
        logging.error(f"Error analyzing gene {gene['gene_id']}: {str(e)}")
        return None

def analyze_variants(genes_df: pd.DataFrame, caspian_variants: Dict[str, List[dict]], 
                    amur_variants: Dict[str, List[dict]], chrom_mapping: Dict[str, str]) -> Tuple[pd.DataFrame, dict]:
    """Analyze all variants and generate summary statistics."""
    logging.info("Starting variant analysis...")
    
    # Create output directory for intermediate results
    os.makedirs('results/intermediate', exist_ok=True)
    
    # Process genes in chunks
    chunk_size = 100  # Smaller chunk size
    results = []
    summary = {
        'total_genes': len(genes_df),
        'genes_with_variants': 0,
        'highly_variable_genes': 0,
        'conserved_genes': 0,
        'total_variation': 0,
        'unique_caspian_total': 0,
        'unique_amur_total': 0,
        'total_snps': 0,
        'total_indels': 0
    }
    
    # Process genes in chunks
    for i in range(0, len(genes_df), chunk_size):
        chunk_genes = genes_df.iloc[i:i+chunk_size]
        chunk_results = []
        
        for _, gene in tqdm(chunk_genes.iterrows(), desc=f"Processing genes {i}-{i+chunk_size}"):
            result = analyze_gene(gene, caspian_variants, amur_variants, chrom_mapping)
            if result:
                chunk_results.append(result)
                
                # Update summary statistics
                if result['total_variants'] > 0:
                    summary['genes_with_variants'] += 1
                if result['variation_per_kb'] > 1:
                    summary['highly_variable_genes'] += 1
                if result['total_variants'] == 0:
                    summary['conserved_genes'] += 1
                
                summary['total_variation'] += result['variation_per_kb']
                summary['unique_caspian_total'] += result['caspian_variants']
                summary['unique_amur_total'] += result['amur_variants']
                summary['total_snps'] += result['caspian_snps'] + result['amur_snps']
                summary['total_indels'] += result['caspian_indels'] + result['amur_indels']
        
        # Save chunk results immediately
        if chunk_results:
            chunk_df = pd.DataFrame(chunk_results)
            chunk_file = f'results/intermediate/chunk_{i}.csv'
            chunk_df.to_csv(chunk_file, index=False)
            logging.info(f"Saved chunk {i} to {chunk_file}")
            results.extend(chunk_results)
        
        # Force garbage collection
        import gc
        gc.collect()
        
        # Log progress
        logging.info(f"Processed {i+len(chunk_genes)}/{len(genes_df)} genes")
    
    # Combine all results
    logging.info("Combining results from all chunks...")
    results_df = pd.concat([pd.read_csv(f) for f in sorted(glob.glob('results/intermediate/chunk_*.csv'))])
    
    # Calculate final summary statistics
    summary['average_variation'] = summary['total_variation'] / len(genes_df)
    summary['median_variation'] = results_df['variation_per_kb'].median()
    summary['max_variation'] = results_df['variation_per_kb'].max()
    summary['shared_variants_total'] = min(summary['unique_caspian_total'], summary['unique_amur_total'])
    
    logging.info("Variant analysis completed successfully")
    return results_df, summary

def plot_variation_distribution(results_df: pd.DataFrame, output_dir: str):
    """Create plots showing variation distribution."""
    logging.info("Generating variation distribution plots...")
    
    try:
        # Set style
        plt.style.use('seaborn')
        
        # Plot 1: Distribution of variation rates
        plt.figure(figsize=(10, 6))
        sns.histplot(data=results_df, x='variation_per_kb', bins=50)
        plt.title('Distribution of Variation Rates (per kb)')
        plt.xlabel('Variation Rate (variants/kb)')
        plt.ylabel('Number of Genes')
        plt.savefig(os.path.join(output_dir, 'variation_distribution.png'))
        plt.close()
        
        # Plot 2: Top variable genes
        top_genes = results_df.nlargest(20, 'variation_per_kb')
        plt.figure(figsize=(12, 8))
        sns.barplot(data=top_genes, x='name', y='variation_per_kb')
        plt.title('Top 20 Most Variable Genes')
        plt.xlabel('Gene Name')
        plt.ylabel('Variation Rate (variants/kb)')
        plt.xticks(rotation=45, ha='right')
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, 'top_variable_genes.png'))
        plt.close()
        
        # Plot 3: Variant type distribution
        variant_types = pd.DataFrame({
            'Subspecies': ['Caspian'] * 2 + ['Amur'] * 2,
            'Type': ['SNP', 'INDEL'] * 2,
            'Count': [
                results_df['caspian_snps'].sum(),
                results_df['caspian_indels'].sum(),
                results_df['amur_snps'].sum(),
                results_df['amur_indels'].sum()
            ]
        })
        
        plt.figure(figsize=(10, 6))
        sns.barplot(data=variant_types, x='Subspecies', y='Count', hue='Type')
        plt.title('Distribution of Variant Types')
        plt.xlabel('Subspecies')
        plt.ylabel('Count')
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, 'variant_types.png'))
        plt.close()
        
        logging.info("Plots generated successfully")
        
    except Exception as e:
        logging.error(f"Error generating plots: {str(e)}")
        raise

def main():
    try:
        # Load chromosome mapping
        chrom_mapping = load_chromosome_mapping('chromosome_mapping.txt')
        
        # Parse GFF file
        genes_df = parse_gff_file('data/reference/panthera_tigris.gff')
        
        # Get variants from VCF files
        caspian_variants = get_variants('data/caspian/aligned/variants.vcf.gz', list(genes_df['chromosome'].unique()))
        amur_variants = get_variants('data/amur/aligned/variants.vcf.gz', list(genes_df['chromosome'].unique()))
        
        # Analyze variants
        results_df, summary = analyze_variants(genes_df, caspian_variants, amur_variants, chrom_mapping)
        
        # Save results
        results_df.to_csv(os.path.join('results', 'tiger_comparison.csv'), index=False)
        
        # Quantitative report (joint VCF + cached pipeline text + gene table below)
        _scripts_dir = os.path.dirname(os.path.abspath(__file__))
        if _scripts_dir not in sys.path:
            sys.path.insert(0, _scripts_dir)
        from analysis_summary_report import write_comprehensive_analysis_summary

        write_comprehensive_analysis_summary(gene_summary=summary)
        
        # Generate plots
        plot_variation_distribution(results_df, 'results')
        
        logging.info("Analysis completed successfully")
        
    except Exception as e:
        logging.error(f"Error in main analysis: {str(e)}")
        raise

if __name__ == "__main__":
    main() 