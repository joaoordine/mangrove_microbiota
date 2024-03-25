# Microbiota analyses 
- Code adapted from Zaramela and Zorz (https://github.com/Microbial-Ecosystems-Lab/asthma_obesity and https://github.com/jkzorz/SituSeq/)

Open tidyverse environment (using conda), then open R in terminal --> install required packages; when finished, open Rstudio, also in terminal; This workflow usually works better than simply oppening Rstudio 

## Go to script 1. Preprocessing 

## Go to script 2. Database creation and Reads Alignment

## Go to script 3. Assigning taxonomy to nanopore reads
## Go to script 3.1 Calculate annotation frequency  

## Go to script 4. Diversity analysis

## Go to script 5. Visualize taxonomy 

## Go to script 7. Correlation heatmap 

Some mangrove isolates were selected for whole-genome sequencing (Illumina NovaSeq 6000). The assembly, annotation and all including analyses can be assessed in the following scripts.
- bioinformatic pipeline adapted from Borelli, et al. "Combining functional genomics and whole-genome sequencing to detect antibiotic resistance genes in bacterial strains co-occurring simultaneously in a Brazilian hospital." Antibiotics 10.4 (2021): 419.

## Go to script 8. Processing raw sequencing files

## Go to script 9. Quality check and filtering

## Go to script 10. Genome assembly and QC

## Go to script 11. Genome annotation

Taxonomical identification of sequenced genomes - TYGS: https://tygs.dsmz.de

Annotation and analysis of secondary metabolite biosynthesis gene clusters was performed with antiSMASH: https://antismash.secondarymetabolites.org/#!/start 

When trying to submit the assembled genomes on Genome - NCBI database, I had to trim the scaffolds to avoid issues. Along with that, I also used an NCBI tool to check for the presence of Illumina adaptors still present in the assembly. For that:
## Go to script 12. NCBI submission 

The bioremediation assessment was performed by analyzing growth curves. The script used to plot growth curves and growth rates of each specific isolate can be seen at
## Go to script 13. Visualize isolates growth
