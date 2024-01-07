# Comparative_genomics_snakes

Files for the comparative analyses of several snake genomes

Starting data include the genomes of three snake species, Sistrurus catenatus, Sistrurus miliarius, and Sistrurus tergeminus. These species genomes have been assembled using a combination of illumina short read data and Pacbio CLR data which have been assembled with MaSuRCA. They were each subsequently annotated with venom gland RNAseq and publicly available protein data using the programs RepeatMasker and MAKER. Annotations were checked and modified manually with mapped RNAseq data.

Following the accumulation of these genomic resources expression was estimated by mapping reads with HiSat2 with expression counts estimated with FeatureCounts. Orthologous groups of genes were inferred with OrthoFinder. Gene trees were inferred for each orthologus group (mafft alignment, iqtree), then reconciled with generax.

This repository contains scripts for these processing steps as well as the subsequent statistical analysis.

Starting from genomes and gffs:

Module 1: Read mapping. Run on a 'per species' basis to facilitate parallelizing across species
Module 2: CDS extraction, orthofinder, creating directory structure, and writing sequence to orthogroup maps for all orthogroups. Run once for all species.
Module 3: Alignment construction, initial tree inference, reconciliation, paralog classification, dissimilarity assessment. Run on a 'per orthogroup' basis to facilitate parallelizing across orthogroups.
Module 4: Summary statistics. Not written yet, but will collect information across orthogroups for some statistical inference.

Part 1: Prep orthofinder. Would be run as a single script for all samples.
1st extract cds - Extract_gff_feature.py
1.5 filter cds
2nd translate cds - translate.py
3 - make orthofinder directory and put faa's there
4 - run orthofinder
5 - read mapping
6 - AA kmer extraction

Part 2: Parse orthofinder output, parallelize across orthogroups.
5 - read orthologs.tsv,  


