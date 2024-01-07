#!/usr/bin/env python

# Additional software necessary to run this:

# (1) biopython
# (2) pandas
# (3) dfply
# (4) mafft
# (5) iqtree
# (6) generax v.XXX

import argparse
from asyncio import subprocess
from ctypes import alignment
import sys, os, shutil
import subprocess as sp
import datetime as dt
import csv
#import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqIO import FastaIO
from dfply import *

########################################
############### ARGUMENTS ##############
########################################

parser = argparse.ArgumentParser(description='')
parser.add_argument("-s","--scripts",
					type=str,
                    default="../Comparative_genomics_snakes.git/scripts/bin",
					help="Directory containing associated scripts")
parser.add_argument("-st","--species_tree",
					type=str,
                    default="0_Inputs/species_tree.nwk",
					help="Species tree file for orthofinder. Newick format")
parser.add_argument("-gd","--genomes_dir",
					type=str,
                    default="0_Inputs/genomes",
					help="Directory containing the genomes for focal taxa. Default is 0_Inputs/genomes")
parser.add_argument("-a","--annotations_dir",
					type=str,
                    default="0_Inputs/annotations",
					help="Directory containing the genome annotatiosn for focal taxa. Default is 0_Inputs/annotations")
parser.add_argument("-otn","--other_taxa_fna",
					type=str,
                    default="0_Inputs/other_taxa_fna",
					help="Directory containing the CDS nucleotide sequences for other taxa")
parser.add_argument("-ota","--other_taxa_faa",
					type=str,
                    default="0_Inputs/other_taxa_faa",
					help="Directory containing the AA sequences for other taxa")
parser.add_argument("-ed","--expression_dir",
					type=str,
                    default="1_read_mapping",
					help="Directory containing readmapping and featureCounts output for each species. Default is '1_read_mapping'")
parser.add_argument("-fna","--fna_dir",
					type=str,
                    default="2.1_CDS_fna",
					help="Directory containing the CDS nucleotide sequences for focal taxa")
parser.add_argument("-faa","--faa_dir",
					type=str,
                    default="2.2_CDS_faa",
					help="Directory containing the AA sequences for focal taxa")
parser.add_argument("-of","--orthofinder_dir",
					type=str,
					default="2.3_Orthofinder",
					help="Directory to write amino acid CDS files for each species. This directory will be written")
parser.add_argument("-md","--map_dir",
					type=str,
                    default="3.1_Species_maps",
					help="Directory containing sequence")
parser.add_argument("-sd","--sequence_dir",
					type=str,
                    default="3.2_OG_sequences",
					help="Directory of where to write a fasta of collected orthogroup sequences")
parser.add_argument("-ad","--alignment_dir",
					type=str,
                    default="3.3_OG_alignments",
					help="Directory of where to write an aligned fasta of sequences")
parser.add_argument("-td","--tree_dir",
					type=str,
                    default="3.4_OG_trees",
					help="Directory of where to run iqtree")
parser.add_argument("-rd","--reconciliations_dir",
					type=str,
                    default="3.5_OG_reconciliations",
					help="Directory of where to run generax to reconcile the genetree")
parser.add_argument("-pcd","--paralog_classifications_dir",
					type=str,
                    default="3.6_OG_classifications",
					help="Directory of where to write the output of orthocaller which generates stats input files")
parser.add_argument("-akd","--AA_kmers_dir",
					type=str,
                    default="3.7_OG_AA_kmers",
					help="Directory of where to write the output of orthocaller which generates stats input files")
parser.add_argument("-oe","--orthogroup_expression_dir",
					type=str,
                    default="3.8_OG_expression",
					help="Directory of where to write the output of expression orthogroup formatter. Will produce an orthogroup specific directory and csv file ascribing expression to genes within an orthogroup based on classifications. Default is '3.8_OG_expression'")
parser.add_argument("-go","--orthogroup_gene_order_dir",
					type=str,
                    default="3.9_OG_geneorder",
					help="Directory of where to write orthogroup gene orders. Will produce an orthogroup specific directory and csv file giving geneorders. Default is '3.9_OG_geneorder'")
parser.add_argument("-dod","--dissimilarity_output_dir",
					type=str,
                    default="4.1_OG_Routputs",
					help="Directory where to write dissimilarity output directory and files. Default is 4.1_OG_Routputs")
parser.add_argument("-cm","--codeml_dir",
					type=str,
                    default="4.2_OG_codeml",
					help="Directory where codeml inference will be conducted. Default is 4.2_OG_codeml")
parser.add_argument("-I","--orthofinder_I",
					type=float,
					default=1.5,
					help="I value to pass to orthofinder. Default is 1.5. I have found 1.3 useful for clustering venom gene families.")
parser.add_argument("-m","--min_focal_taxa",
					type=int,
                    default=2,
					help="Minimum number of focal taxa to process an orthogroup downstream. Default is 2")
parser.add_argument("-ct","--counts_tpm",
					type=str,
                    default="tpm",
					help="string that identifies whether you would like your ouput as read counts or transcripts per million (tpm). Options are 'count' or 'tpm'. TPM is default")
parser.add_argument("-iq","--iqtree",
					type=str,
                    default="iqtree",
					help="Path to iqtree. Default assumes it is in path")
parser.add_argument("-t","--threads",
					type=int,
					default=1,
					help="Number of processing threads. (default: %(default)s)")
args = parser.parse_args()


########################################
################# SETUP ################
########################################

genomes_dir = args.genomes_dir
annotations_dir = args.annotations_dir
map_dir = args.map_dir
fna_dir = args.fna_dir
faa_dir = args.faa_dir
orthofinder_dir = args.orthofinder_dir
otn_dir = args.other_taxa_fna
ota_dir = args.other_taxa_faa
sequence_dir = args.sequence_dir
alignment_dir = args.alignment_dir
tree_dir = args.tree_dir
reconciliations_dir= args.reconciliations_dir
paralog_classifications_dir = args.paralog_classifications_dir
AA_kmers_dir = args.AA_kmers_dir
orthogroup_gene_order_dir = args.orthogroup_gene_order_dir
dissimilarity_output_dir = args.dissimilarity_output_dir
codeml_dir = args.codeml_dir
scripts = args.scripts
species_tree = args.species_tree
expression_dir = args.expression_dir
orthogroup_expression_dir = args.orthogroup_expression_dir
orthofinder_I = args.orthofinder_I 
min_focal_taxa=args.min_focal_taxa
counts_tpm = args.counts_tpm
iqtree = args.iqtree
threads = args.threads

########################################
############## FUNCTIONS ###############
########################################

def run_make_directory_structure(fna_dir, faa_dir, orthofinder_dir, map_dir, sequence_dir, alignment_dir, tree_dir, reconciliations_dir, paralog_classifications_dir, AA_kmers_dir, orthogroup_expression_dir, dissimilarity_output_dir):
    mkdir_call = "mkdir " + fna_dir + " " + faa_dir + " 2.2.2_CDS_stops " + orthofinder_dir + " " + map_dir + " " + sequence_dir + " " + alignment_dir + " " + tree_dir + " " + reconciliations_dir + " " + paralog_classifications_dir + " " + AA_kmers_dir + " " + orthogroup_expression_dir + " " + dissimilarity_output_dir + " " + orthogroup_gene_order_dir + " " + codeml_dir
    sp.call(mkdir_call, shell=True)


def run_prep_orthofinder(scripts, genomes_dir, annotations_dir, fna_dir, faa_dir, expression_dir, orthofinder_dir, orthofinder_I, threads):
    prep_of_call = "python " + scripts + "/prep_orthofinder.py -g " + genomes_dir + " -a " + annotations_dir + " -s " + scripts + " -cd " + fna_dir + " -aa " + faa_dir + " -rm " + expression_dir + " -of " + orthofinder_dir + " -I " + str(orthofinder_I) + " -t " + str(threads)
    sp.call(prep_of_call, shell=True)


def run_process_orthogroups(orthofinder_dir, scripts, genomes_dir, min_focal_taxa, map_dir):
    orthofinder_results_folder =  orthofinder_dir + "/OrthoFinder"
    orthofinder_input = os.listdir(orthofinder_results_folder)[0]
    orthofinder_results_input = orthofinder_results_folder + "/" + orthofinder_input
    process_orthogroups_call = "python " + scripts + "/process_orthogroups.py -of " + orthofinder_results_input + " -g " + genomes_dir + " -m " + str(min_focal_taxa) + " -out " + map_dir
    sp.call(process_orthogroups_call, shell=True)


########################################
################# CODE #################
########################################
print(scripts)
run_make_directory_structure(fna_dir, faa_dir, orthofinder_dir, map_dir, sequence_dir, alignment_dir, tree_dir, reconciliations_dir, paralog_classifications_dir, AA_kmers_dir, orthogroup_expression_dir, dissimilarity_output_dir)

print(1)
run_prep_orthofinder(scripts, genomes_dir, annotations_dir, fna_dir, faa_dir, expression_dir, orthofinder_dir, orthofinder_I, threads)
print(2)
run_process_orthogroups(orthofinder_dir, scripts, genomes_dir, min_focal_taxa, map_dir)

