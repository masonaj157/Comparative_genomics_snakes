#!/usr/bin/env python

# Additional software necessary to run this:

# (1) HiSat2
# (2) featureCount

import argparse
from asyncio import subprocess
import sys, os, shutil
import subprocess as sp
import datetime as dt
import csv
import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqIO import FastaIO
from dfply import *

########################################
############### ARGUMENTS ##############
########################################

parser = argparse.ArgumentParser(description='')
parser.add_argument("-s","--species",
					type=str,
					help="Name of focal species. Should match genome file, annotations file, and a species assignment listed in the sample mapping file")
parser.add_argument("-g","--genome_dir",
					type=str,
					help="Directory containing the input genomes as fastas")
parser.add_argument("-a","--annotations_dir",
					type=str,
					help="Directory containing the input annotations as gff3")
parser.add_argument("-r","--reads_dir",
					type=str,
					help="Directory containing the input fastqs from RNAseq data each named by sample. Forward and reverse reads should be designated with \'_F\' and \'_R\' as in (Sample_F.fastq).")
parser.add_argument("-m","--species_map",
					type=str,
					help="Comma separated file file where each lines lists a species and its species assignment with headers \'Sample,Species\'")
parser.add_argument("-hp","--hisat2_path",
					type=str,
                    default='',
					help="Directory path to HiSat2 and associated python scripts. Default assumes it is in your path")
parser.add_argument("-p","--processes",
					type=int,
                    default=20,
					help="Number of processes to be run by HiSat2. Default is 20")
parser.add_argument("-o","--output_directory",
					type=str,
					help="Directory of where to write a \'Species\' directory and the subsequent read mapping files.")
args = parser.parse_args()

########################################
################# SETUP ################
########################################

species = args.species
genome_dir = args.genome_dir
annotations_dir = args.annotations_dir
reads_dir = args.reads_dir
species_map = args.species_map
hisat2_path = args.hisat2_path
output_dir = args.output_directory
procs = args.processes


if hisat2_path != '' and hisat2_path[-1] != '/':
    hisat2_path = hisat2_path + '/'

print("""
Running HiSat2, processing output alignment files, and estimating per transcript counts with featureCounts
""")

print("\n"+dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" ::: Starting...")
print("\tSpecies --> "+ species)
print("\tGenome directory --> "+ genome_dir)
print("\tAnnotations directory --> "+ annotations_dir)
print("\tReads directory --> "+ reads_dir)
print("\tSample to species mapping file --> "+ species_map)
print("\tOutput directory --> "+ output_dir)

########################################
############## FUNCTIONS ###############
########################################

def make_samples_list(species,species_map):
    samples_species = pd.read_csv(species_map, ',')
    relevant_samples = samples_species[samples_species["Species"] == species]
    samples = list(relevant_samples["Sample"])
    return(samples)

def genome_preprocessing(species, genome_dir, annotations_dir, output_dir, hisat2_path):
    convert_gff = "gffread " + annotations_dir + '/' + species + ".gff -T -o " + output_dir + "/" + species + ".gtf"
    sp.call(convert_gff, shell=True)
    #print(convert_gff)
    if hisat2_path != '':
        #print(hisat2_path)
        generate_ss = "python "+ hisat2_path + "hisat2_extract_splice_sites.py " + output_dir + "/" + species + ".gtf > "+ output_dir + "/" + species + ".ss"
        sp.call(generate_ss, shell=True)
        #print(generate_ss)
        generate_exon = "python "+ hisat2_path + "hisat2_extract_exons.py " + output_dir + "/" + species + ".gtf > "+ output_dir + "/" + species + ".exon"
        sp.call(generate_exon, shell=True)
        #print(generate_exon)
        build_index = hisat2_path + "hisat2-build -p " + str(procs) + " --ss " + output_dir + "/" + species + ".ss --exon " + output_dir + "/" + species + ".exon -f " + genome_dir + "/" + species + ".fasta " + output_dir + "/" + species + "_annotations"
        sp.call(build_index, shell=True)
        #print(build_index)
    else:
        #print(hisat2_path)
        generate_ss = "python hisat2_extract_splice_sites.py " + output_dir + "/" + species + ".gtf > "+ output_dir + "/" + species + ".ss"
        sp.call(generate_ss, shell=True)
        #print(generate_ss)
        generate_exon = "python hisat2_extract_exons.py " + output_dir + "/" + species + ".gtf > "+ output_dir + "/" + species + ".exon"
        sp.call(generate_exon, shell=True)
        #print(generate_exon)
        build_index = "hisat2-build -p " + str(procs) + " --ss " + output_dir + "/" + species + ".ss --exon " + output_dir + "/" + species + ".exon -f " + genome_dir + "/" + species + ".fasta " + output_dir + "/" + species + "_annotations"
        sp.call(build_index, shell=True)
        #print(build_index)

def run_Hisat(samples,species,reads_dir,output_dir):
    for sample in samples:
        mkdir_command = "mkdir " + output_dir + "/" + species + '/' + sample 
        #print(mkdir_command)
        sp.call(mkdir_command, shell=True)
        hisat = hisat2_path + "hisat2 -p " + str(procs) + " -k 10 --dta -x " + output_dir + "/" + species + "_annotations -1 " + reads_dir + "/" + sample + "_F.fastq.gz -2 " + reads_dir + "/" + sample + "_R.fastq.gz -S " + output_dir + "/" + species + "/" + sample + "/" + sample + "_aln.sam"
        sp.call(hisat, shell=True)
        #print(hisat)

def convert_to_bam(samples, species, output_dir):
    for sample in samples:
        samtools_call = "samtools sort -@ 16 -o " + output_dir + "/" + species + "/" + sample + "/" + sample + "_aln.bam " + output_dir + "/" + species + "/" + sample + "/" + sample + "_aln.sam"
        #print(samtools_call)
        sp.call(samtools_call,shell=True)
        samtools_index = "samtools index " + output_dir + "/" + species + "/" + sample + "/" + sample + "_aln.bam"
        #print(samtools_index)
        sp.call(samtools_index,shell=True)


def run_featureCounts(samples, species, output_dir, procs):
    samples_string = ''
    for sample in samples:
        samples_string = samples_string + " " + output_dir + "/" + species + "/" + sample + "/" + sample + "_aln.bam"
    fc_command = "featureCounts -a " + output_dir + "/" + species + ".gtf -o " + output_dir + "/" + species + "/" + species + "_counts -T " + str(procs)+ " -p -C -g transcript_id -t CDS" + samples_string
    #print(fc_command)
    sp.call(fc_command, shell=True)

########################################
################# CODE #################
########################################

samples = make_samples_list(species, species_map)

mk_species_dir = "mkdir " + output_dir + "/" + species
print(mk_species_dir)
sp.call(mk_species_dir, shell=True)

genome_preprocessing(species, genome_dir, annotations_dir, output_dir,hisat2_path)

run_Hisat(samples,species,reads_dir,output_dir)

convert_to_bam(samples, species, output_dir)

run_featureCounts(samples, species, output_dir, procs)