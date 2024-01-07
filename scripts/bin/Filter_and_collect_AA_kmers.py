#!/usr/bin/env python

import argparse
import csv
import pandas as pd
import subprocess as sp
import sys, os, shutil
from Bio import SeqIO

# Command line options
parser = argparse.ArgumentParser(description='This script will take an orthogroup, directory of amino acid fastas, and a sequence map and for each species 1) read the fasta 2) filter to relevant sequences like orthogroups 3) calculate kmers and 4) output the table as a csv file.')
parser.add_argument("-og","--orthogroup",
					type=str,
					help="Directory containing the amino acid fastas for each species. Files named as 'Species.faa'")
parser.add_argument("-gd","--genomes_dir",
					type=str,
                    default="0_Inputs/genomes",
					help="Directory containing the genomes for each species. Files named as 'Species.fasta'")
parser.add_argument("-ad","--amino_acid_dir",
					type=str,
                    default="2.2_CDS_faa",
					help="Directory containing the amino acid fastas for each species. Files named as 'Species.faa'")
parser.add_argument("-md","--map_dir",
					type=str,
                    default="3.1_Species_maps",
					help="Directory containing the orthogroup map'")
parser.add_argument("-k","--kmer_size",
					type=int,
                    default=20,
					help="AA kmer size for counting")                    
parser.add_argument("-o","--output_dir",
                    type=str,
                    default="3.7_OG_AA_kmers",
                    help="name of output directory of where to write the folders and output")
args=parser.parse_args()

orthogroup=args.orthogroup
genomes_dir=args.genomes_dir
amino_acid_dir=args.amino_acid_dir
map_dir=args.map_dir
kmer_size=args.kmer_size
output_dir=args.output_dir

#########################################################################################

def find_species(genomes_dir) :
	species_list = list(os.listdir(genomes_dir))
	species_list = [x.split('.')[0] for x in species_list]
	return(species_list)


def extract_kmers(seq, kmer_size):
    if len(seq.seq) <= kmer_size:
        print(seq.id)
        kmers = 'NA'
    else:
        kmers = []
        start = 0
        end = start + kmer_size
        while end <= len(seq.seq):
            kmers.append(seq.seq[start:end])
            start +=1
            end = start + kmer_size
    return(kmers)


def extract_kmers(seq, kmer_size):
    if len(seq.seq) <= kmer_size:
        print(seq.id)
        kmers = [seq.id,'NA']
    else:
        kmers = []
        start = 0
        end = start + kmer_size
        while end <= len(seq.seq):
            kmers.append([seq.id,seq.seq[start:end]])
            start +=1
            end = start + kmer_size
    return(kmers)


def collect_all_kmers(sequences, kmer_size):
    all_kmers = []
    for seq in sequences:
        kmers = extract_kmers(seq, kmer_size)
        for kmer in kmers:
            all_kmers.append(kmer)
    
    return(all_kmers)



def filter_seqs(sequences,map,species):
    relevant_seq_names = list(map[map["Species"] == species]["Sequence"])
    #print(relevant_seq_names)
    relevant_seqs = [seq for seq in sequences if seq.id in relevant_seq_names]
    return(relevant_seqs)


#########################################################################################

species_list = find_species(genomes_dir)
if '' in species_list:
    species_list.remove('')
if '.dir' in species_list:
    species_list.remove('.dir')

OG_output_dir = "mkdir " + output_dir + "/" + orthogroup
sp.call(OG_output_dir, shell=True)

map_file = map_dir + "/" + orthogroup + '.map'
map = pd.read_csv(map_file,",")

for species in species_list:
    amino_acid_fasta = amino_acid_dir + "/" + species + ".faa"
    sequences = list(SeqIO.parse(amino_acid_fasta,"fasta"))

    sequences = filter_seqs(sequences,map,species)

    all_kmers = collect_all_kmers(sequences, kmer_size)
    

    output_file = output_dir + "/" + orthogroup + "/" + species + "_" + orthogroup + "_AA_kmers.csv"

    with open(output_file,'w') as outFile:
        outfile_writer = csv.writer(outFile, delimiter = ',')
        outfile_writer.writerow(['seqid','kmer'])
        for entry in all_kmers:
            outfile_writer.writerow(entry)

    outFile.close()

