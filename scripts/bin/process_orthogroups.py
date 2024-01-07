#!/usr/bin/env python

# Additional software necessary to run this:

# (1) biopython
# (2) pandas
# (3) dfply
# (4) mafft
# (5) iqtree
# (6) generax v.XXX

import argparse
from ast import keyword
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
#parser.add_argument("-og","--orthogroup",
#					type=str,
#					help="Orthogroup from an OrthoFinder orthogroups.tsv file")
parser.add_argument("-of","--orthofinder_input",
					type=str,
					help="OrthoFinder results folder containing the output of an orthofinder run")
parser.add_argument("-g","--genomes_dir",
					type=str,
					help="Directory containing input genomes. Used to differentiate focal species from non-focal taxa")
parser.add_argument("-m","--min_focal_taxa",
					type=int,
                    default=2,
					help="Minimum number of focal taxa to process an orthogroup downstream. Default is 2")
parser.add_argument("-","--keywords_file",
					type=str,
                    default='0_Inputs/keywords.txt',
					help="Text file containing keywords to find among orthogroups")
parser.add_argument("-out","--output_directory",
					type=str,
                    default = "3.1_Species_maps",
					help="Directory of where to write final map files. Default is 3.1_Species_maps")
parser.add_argument("-ko","--keywords_output",
					type=str,
                    default = "0_Inputs/keyword_orthogroups.csv",
					help="Directory of where to write final map files. Default is 3.1_Species_maps")
args = parser.parse_args()

########################################
################# SETUP ################
########################################

#orthogroup = args.orthogroup
orthofinder_input = args.orthofinder_input
genomes_dir = args.genomes_dir
out = args.output_directory
min_focal_taxa = args.min_focal_taxa
keywords_file = args.keywords_file
keywords_output = args.keywords_output

#input = os.path.abspath(args.input)
#input_name = os.path.basename(input)
#input_name2 = input_name.split(".fasta")[0]
#reads = os.path.abspath(args.reads)
#cov = args.cov
#prop = args.prop

print("""
Processing output of OrthoFinder.
""")

print("\n"+dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" ::: Starting...")
#print("\tOrthogroup --> "+ orthogroup)
print("\tOutput directory --> "+ out)

########################################
############## FUNCTIONS ###############
########################################

def find_species(genome_dir) :
    species_list = list(os.listdir(genome_dir))
    species_list = [x.split('.fasta')[0] for x in species_list]
    if ".DS_Store" in species_list:
        species_list.remove(".DS_Store")
    return(species_list)


def test_for_focal_taxa(tsv_df, species_list, min_focal_taxa):
    lose_list = []
    for og in list(tsv_df["Orthogroup"]):
        og_row = tsv_df[tsv_df["Orthogroup"] == og]
        counter = 0
        for species in species_list:
            test_val = og_row[species]
            if not test_val.isnull().values.any():
                counter += 1
        if counter < min_focal_taxa:
            lose_list.append(og)
    tsv_df = tsv_df[~tsv_df["Orthogroup"].isin(lose_list)]
    #print(tsv_df["Orthogroup"])
    return(tsv_df)
    


def make_map(og_row):
    seq_map = []
    for species in list(og_row.columns)[1:]:
        species_seqs = og_row[species]
        if not species_seqs.isnull().values.any():
            for seqlist in list(species_seqs.str.split(", ")):
                for seq in seqlist:
                    if " " in seq:
                        seq = seq.split(" ")[0]
                    new_name = species + '_' + seq
                    row = [seq,species,new_name]
                    seq_map.append(row)
    return(seq_map)


def test_for_keywords(tsv_df,keywords):
    keywords_orthogroups = []
    for og in list(tsv_df["Orthogroup"]) :
        test = False
        og_row = tsv_df[tsv_df["Orthogroup"] == og]
        for species in list(og_row.columns)[1:]:
            species_seqs = og_row[species]
            if not species_seqs.isnull().values.any():
                for seqlist in list(species_seqs.str.split(", ")):
                    for seq in seqlist:
                        if " " in seq:
                            seq = seq.split(" ")[0]
                            for key in keywords:
                                if key in seq:
                                    test = True
                                    keywords_orthogroups.append([key,og])
                                    break
                        if test == True:
                            break
                    if test == True:
                        break
            if test == True:
                break
    return(keywords_orthogroups)

    

########################################
################# CODE #################
########################################

species_list = find_species(genomes_dir)

print(species_list)

keywords = []
with open(keywords_file) as kw:
    keywords = kw.readlines()


keywords = [x.split('\n')[0] for x in keywords]
keywords.remove('')


tsv = orthofinder_input + "/Orthogroups/Orthogroups.tsv"

tsv_df = pd.read_csv(tsv,"\t")

tsv_df = test_for_focal_taxa(tsv_df, species_list, min_focal_taxa)

for og in list(tsv_df["Orthogroup"]) :
    og_row = tsv_df[tsv_df["Orthogroup"] == og]
    seq_map = make_map(og_row)
    
    outfile = out + '/' + og + '.map'
    with open(outfile, 'w') as csv_file:
        csv_writer = csv.writer(csv_file, delimiter = ',')
        csv_writer.writerow(['Sequence','Species','New_name'])
        for row in seq_map:
            csv_writer.writerow(row)

    csv_file.close()


keywords_output_tab = test_for_keywords(tsv_df,keywords)

with open(keywords_output, 'w') as csv_file:
    csv_writer = csv.writer(csv_file, delimiter = ',')
    csv_writer.writerow(['keyword','Orthogroup'])
    for row in keywords_output_tab:
        csv_writer.writerow(row)

csv_file.close()
