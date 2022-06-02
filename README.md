# Comparative_genomics_snakes

Files for the comparative analyses of several snake genomes

Starting data include the genomes of three snake species, Sistrurus catenatus, Sistrurus miliarius, and Sistrurus tergeminus. These species genomes have been assembled using a combination of illumina short read data and Pacbio CLR data which have been assembled with MaSuRCA. They were each subsequently annotated with venom gland RNAseq and publicly available protein data using the programs RepeatMasker and MAKER. Annotations were checked and modified manually with mapped RNAseq data.

Following the accumulation of these genomic resources expression was estimated by mapping reads with HiSat2 with expression counts estimated with FeatureCounts. Orthologous groups of genes were inferred with OrthoFinder. Gene trees were inferred for each orthologus group (mafft alignment, iqtree), then reconciled with generax.

This repository contains scripts for these processing steps as well as the subsequent statistical analysis. 


