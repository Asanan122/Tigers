#!/bin/bash

# Reference genome
REF="data/reference/panthera_tigris.fasta"

# Caspian Tiger
echo "Calling variants for Caspian Tiger..."
bcftools mpileup -Ou -f $REF data/caspian/aligned/SRR18572400.sorted.bam | \
    bcftools call -mv -Ov -o data/caspian/aligned/variants.vcf

# Amur Tiger
echo "Calling variants for Amur Tiger..."
bcftools mpileup -Ou -f $REF data/amur/aligned/SRR31485304.sorted.bam | \
    bcftools call -mv -Ov -o data/amur/aligned/variants.vcf

echo "Variant calling complete!" 