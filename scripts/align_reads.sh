#!/bin/bash

# Exit on error
set -e

# Configuration
REFERENCE="data/reference/panthera_tigris.fasta"
THREADS=8

# Create output directories
mkdir -p data/amur/aligned
mkdir -p data/caspian/aligned

# Function to align paired-end reads
align_reads() {
    local sample=$1
    local r1=$2
    local r2=$3
    local output_prefix=$4

    echo "Aligning reads for ${sample}..."
    
    # Align reads using BWA-MEM
    bwa mem -t ${THREADS} \
        ${REFERENCE} \
        ${r1} ${r2} \
        | samtools view -bS - \
        | samtools sort -@ ${THREADS} -o ${output_prefix}.sorted.bam -

    # Index the BAM file
    samtools index ${output_prefix}.sorted.bam
    
    echo "Alignment completed for ${sample}"
}

# Index reference genome for BWA if not already done
if [ ! -f "${REFERENCE}.bwt" ]; then
    echo "Indexing reference genome for BWA..."
    bwa index ${REFERENCE}
fi

# Align Caspian tiger reads
align_reads \
    "Caspian" \
    "data/caspian/SRR18572400_1.fastq" \
    "data/caspian/SRR18572400_2.fastq" \
    "data/caspian/aligned/SRR18572400"

# Align Amur tiger reads
align_reads \
    "Amur" \
    "data/amur/SRR31485304_1.fastq" \
    "data/amur/SRR31485304_2.fastq" \
    "data/amur/aligned/SRR31485304"

echo "All alignments completed!" 