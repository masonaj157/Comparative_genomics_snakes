#!/usr/bin/env python

# Additional software necessary to run this:
# (1) bwa 
# (2) samtools
# (3) bedtools
# (4) biopython
# (5) pandas
# (6) dfply

import argparse
import sys, os, shutil
import subprocess as sp
import datetime as dt
import csv
#import numpy as np
#import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqIO import FastaIO
from dfply import *

########################################
############### ARGUMENTS ##############
########################################

parser = argparse.ArgumentParser(description='')
parser.add_argument("-g","--genome_dir",
					type=str,
					help="Directory containing the input genomes as fastas")
parser.add_argument("-ot","--other_taxa_dir",
					type=str,
					default = "0_Inputs/other_taxa_fna",
					help="Directory containing the .fna files for additional reference taxa")
parser.add_argument("-a","--annotations_dir",
					type=str,
					help="Directory containing the input annotations as gff3")
parser.add_argument("-s","--scripts_dir",
					type=str,
					default='../Comparative_genomics_snakes.git/scripts/bin',
					help="Directory containing processing scripts")
parser.add_argument("-cd","--cds_dir",
					type=str,
					default="2.1_CDS_fna",
					help="Directory to write nucleotide CDS files for each species. This directory will be written")
parser.add_argument("-aa","--aa_dir",
					type=str,
					default="2.2_CDS_faa",
					help="Directory to write amino acid CDS files for each species. This directory will be written")
parser.add_argument("-rm","--readmapping_dir",
					type=str,
					default="1_read_mapping",
					help="Directory where read mapping has been done for each species. Scripts assumes it is organized a dir/Species/Species_counts")
parser.add_argument("-of","--orthofinder_dir",
					type=str,
					default="2.3_Orthofinder",
					help="Directory to write amino acid CDS files for each species. This directory will be written")
parser.add_argument("-I","--orthofinder_I",
					type=float,
					default=1.5,
					help="I value to pass to orthofinder. Default is 1.5. I have found 1.3 useful for clustering venom gene families.")
parser.add_argument("-t","--threads",
					type=int,
					default=1,
					help="Number of processing threads. (default: %(default)s)")
args = parser.parse_args()

########################################
################# SETUP ################
########################################

genome_dir = args.genome_dir
other_taxa_dir = args.other_taxa_dir
annotations_dir = args.annotations_dir
scripts_dir = args.scripts_dir
cds_dir = args.cds_dir
aa_dir = args.aa_dir
orthofinder_dir = args.orthofinder_dir
readmapping = args.readmapping_dir
orthofinder_I = args.orthofinder_I
threads = args.threads

print("""
Transcript Presence/Absence Checker
""")

#print("\n"+dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" ::: Starting...")
#print("\tInput --> "+ input_name)
#print("\tReads --> "+ os.path.basename(reads))
#print("\tThreads --> " + str(threads))

########################################
############## FUNCTIONS ###############
########################################

def find_species(genome_dir) :
	species_list = list(os.listdir(genome_dir))
	if ".DS_Store" in species_list:
		species_list.remove(".DS_Store")
	species_list = [x.split('.fasta')[0] for x in species_list]
	species_list = [x.split('.fna')[0] for x in species_list]
	return(species_list)


def Extract_gff(species_list,genome_dir,annotations_dir) :
	mk_CDS = "mkdir " + cds_dir
	sp.call(mk_CDS, shell=True)
	for species in species_list :
		print(species)
		extract_call = "python " + scripts_dir + "/Extract_gff_feature_v0.2.py -s " + genome_dir + "/" + species + ".fasta -g " + annotations_dir + "/" + species + ".gff -f CDS -o " + cds_dir + "/" + species + ".fna"
		sp.call(extract_call, shell=True)
		string = species + ".fna written to " + cds_dir
		print(string)
	ls_call = "ls " + cds_dir + "/"
	sp.call(ls_call, shell=True)


def Test_CDS_expression_names(species_list) :
	for species in species_list:
		print(species)
		run_test = "python " + scripts_dir + "/Test_expression_CDS_match.py -fc " + readmapping + "/" + species + "/" + species + "_counts -cds " + cds_dir + "/" + species + ".fna -g " + annotations_dir + "/" + species + ".gff"
		sp.call(run_test, shell=True)
	ls_call = "ls " + cds_dir + "/"
	sp.call(ls_call, shell=True)


def Translate(species_list, additional_species):
	mk_CDS_aa = "mkdir " + aa_dir
	sp.call(mk_CDS_aa, shell=True)
	mycwd = os.getcwd()
	os.chdir(cds_dir)
	for species in species_list:
		translate_call = "python ../" + scripts_dir + "/Translate_n_filter.py -f " + species + ".fna -o " + species
		sp.call(translate_call,shell=True)
		move_files_call = "mv " + species + ".faa ../" + aa_dir + "/"
		sp.call(move_files_call,shell=True)
		move_files_call = "mv " + species + "_internal_stops.fasta ../2.2.2_CDS_stops/ "
		sp.call(move_files_call,shell=True)
		move_files_call = "mv " + species + "_no_stop.fasta ../_CDS_stops/ "
		sp.call(move_files_call,shell=True)
		move_files_call = "mv " + species + "_cds_filtered.fasta " + species + ".fna"
		sp.call(move_files_call,shell=True)
	os.chdir(mycwd)
	ls_call = "ls " + cds_dir + "/"
	sp.call(ls_call, shell=True)
	os.chdir(other_taxa_dir)
	mkdir_call = "mkdir ../other_taxa_faa"
	for species in additional_species:
		translate_call = "python ../../" + scripts_dir + "/Translate_n_filter.py -f " + species + ".fna -o " + species
		sp.call(translate_call,shell=True)
		move_files_call = "mv " + species + ".faa ../other_taxa_faa/"
		sp.call(move_files_call,shell=True)
		remove_files_call = "rm " + species + "_internal_stops.fasta"
		sp.call(remove_files_call,shell=True)
		remove_files_call = "mv " + species + "_no_stop.fasta"
		sp.call(remove_files_call,shell=True)
		move_files_call = "mv " + species + "_cds_filtered.fasta " + species + ".fna"
		sp.call(move_files_call,shell=True)
	os.chdir(mycwd)
	ls_call = "ls " + cds_dir + "/"
	sp.call(ls_call, shell=True)



def run_Orthofinder(orthofinder_dir, aa_dir, orthofinder_I, threads):
	mk_ortho = "mkdir " + orthofinder_dir
	sp.call(mk_ortho,shell=True)
	copy_call = "cp " + aa_dir + "/*.faa " + orthofinder_dir + "/"
	sp.call(copy_call,shell = True)
	copy_call = "cp 0_Inputs/other_taxa_faa/*.faa " + orthofinder_dir + "/"
	sp.call(copy_call,shell = True)
	orthofinder_call = "orthofinder -f " + orthofinder_dir + " -og -I " + str(orthofinder_I) + " -t " + str(threads)
	sp.call(orthofinder_call,shell = True)  


########################################
################# CODE #################
########################################

## Create CDS files from genomes and annotations
species_list = find_species(genome_dir)
print(species_list)
additional_species = find_species(other_taxa_dir)

print(1)
Extract_gff(species_list,genome_dir,annotations_dir)
print(2)
Test_CDS_expression_names(species_list)
print(3)
Translate(species_list, additional_species)
print(4)
run_Orthofinder(orthofinder_dir, aa_dir, orthofinder_I, threads)
print(5)
