#!/bin/bash

# Exit on error
set -e

# Configuration
REFERENCE="data/reference/panthera_tigris.fasta"
GATK="./gatk-4.4.0.0/gatk"

# Create directories if they don't exist
mkdir -p data/reference/temp

# Create sequence dictionary for reference if it doesn't exist
if [ ! -f "${REFERENCE}.dict" ]; then
    echo "Creating sequence dictionary..."
    ${GATK} CreateSequenceDictionary \
        -R ${REFERENCE} \
        -O ${REFERENCE}.dict
fi

# Index the reference if it doesn't exist
if [ ! -f "${REFERENCE}.fai" ]; then
    echo "Indexing reference genome..."
    samtools faidx ${REFERENCE}
fi

echo "Reference genome preparation completed!" 