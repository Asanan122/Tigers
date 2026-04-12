#!/bin/bash

# Exit on error
set -e

# Configuration
REFERENCE="data/reference/panthera_tigris.fasta"
GATK="./gatk-4.4.0.0/gatk"
THREADS=8
MEMORY="32G"

# Create output directories
mkdir -p data/vcf
mkdir -p data/analysis

# Function to process each sample
process_sample() {
    local sample=$1
    local bam_file="data/${sample}/${sample}.sorted.bam"
    local output_prefix="data/vcf/${sample}"
    
    echo "Processing sample: ${sample}"
    
    # Mark duplicates
    ${GATK} MarkDuplicates \
        -I ${bam_file} \
        -O ${output_prefix}.dedup.bam \
        -M ${output_prefix}.metrics.txt \
        --CREATE_INDEX true \
        --VALIDATION_STRINGENCY SILENT
    
    # Create BAM index
    samtools index ${output_prefix}.dedup.bam
    
    # Create GVCF for joint calling
    ${GATK} HaplotypeCaller \
        -R ${REFERENCE} \
        -I ${output_prefix}.dedup.bam \
        -O ${output_prefix}.g.vcf.gz \
        -ERC GVCF \
        --native-pair-hmm-threads ${THREADS} \
        --tmp-dir /tmp
}

# Process each sample
process_sample "caspian/SRR18572400"
process_sample "amur/SRR31485304"

# Joint genotyping
${GATK} GenotypeGVCFs \
    -R ${REFERENCE} \
    -V data/vcf/caspian/SRR18572400.g.vcf.gz \
    -V data/vcf/amur/SRR31485304.g.vcf.gz \
    -O data/vcf/joint_called.vcf.gz \
    --tmp-dir /tmp

# Apply basic filtering
${GATK} VariantFiltration \
    -R ${REFERENCE} \
    -V data/vcf/joint_called.vcf.gz \
    -O data/vcf/joint_called.filtered.vcf.gz \
    --filter-expression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" \
    --filter-name "basic_filters" \
    --tmp-dir /tmp

# Generate basic statistics
${GATK} VariantsToTable \
    -R ${REFERENCE} \
    -V data/vcf/joint_called.filtered.vcf.gz \
    -F CHROM -F POS -F REF -F ALT -F QUAL -F FILTER -F AC -F AF -F AN -F BaseQRankSum -F ClippingRankSum -F DP -F FS -F MQ -F MQRankSum -F QD -F ReadPosRankSum -F SOR \
    -O data/analysis/variant_stats.txt \
    --tmp-dir /tmp

echo "Variant calling completed for all samples" 