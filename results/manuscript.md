# Genetic Analysis of Caspian and Amur Tiger Subspecies

## Abstract
This study presents a comprehensive genetic analysis comparing the Caspian (Panthera tigris virgata) and Amur (Panthera tigris altaica) tiger subspecies. Using whole-genome sequencing data, we identified and characterized genetic variants to understand the evolutionary relationships and unique adaptations of these subspecies.

## Methods
### Data Processing Pipeline
1. Quality Control: Raw sequencing reads were processed using FastQC and trimmed for quality
2. Alignment: Reads were aligned to the Panthera tigris reference genome using BWA
3. Variant Calling: GATK HaplotypeCaller was used to identify genetic variants
4. Analysis: Statistical analysis and visualization of variants was performed using Python

## Results

### Variant Calling Results
The variant calling pipeline identified a total of 2,169,502 genetic variants across both tiger subspecies. The distribution of variants is as follows:
- Total Variants: 2,169,502
- Shared Variants: 73,119 (3.37%)
- Caspian Unique Variants: 2,008,535 (92.59%)
- Amur Unique Variants: 87,848 (4.04%)
- Homozygous Differences: 194 (0.009%)

### Statistical Analysis
The analysis revealed significant genetic differentiation between the two subspecies:
1. The Caspian tiger shows a much higher number of unique variants (2,008,535) compared to the Amur tiger (87,848)
2. Only 3.37% of variants are shared between the subspecies
3. There are 194 positions where the tigers are fixed for different alleles (homozygous differences)

### Visualization and Interpretation
Two key visualizations were generated to understand the genetic differences:

1. Variant Density Plot (`variant_density.png`)
   - Shows the distribution of variants across chromosomes
   - Helps identify regions of high genetic divergence
   - Useful for identifying potential regions under selection

2. Variant Distribution Plot (`variant_distribution.png`)
   - Illustrates the types and frequencies of genetic variants
   - Helps understand the nature of genetic differences
   - Useful for identifying patterns of evolution

## Discussion
The large number of unique variants in the Caspian tiger sample compared to the Amur tiger sample suggests significant genetic divergence between these subspecies. This could be due to:
1. Geographic isolation and different evolutionary pressures
2. Different adaptation strategies to their respective environments
3. Potential differences in population history and effective population sizes

The low number of shared variants (3.37%) indicates that these subspecies have been genetically isolated for a significant period, supporting their classification as distinct subspecies.

## Conclusions
This study provides valuable insights into the genetic differences between Caspian and Amur tigers. The findings support the distinct nature of these subspecies and highlight the importance of conservation efforts for both populations. Further research is needed to:
1. Understand the functional significance of the identified variants
2. Investigate the evolutionary history of these subspecies
3. Develop conservation strategies based on genetic diversity

## References
1. GATK Best Practices for Variant Calling
2. BWA Documentation
3. FastQC Documentation 