#!/usr/bin/env python

import argparse
import csv
from Bio import SeqIO

parser = argparse.ArgumentParser(description='Takes a fasta alignment and mapping file of individuals\' species assignments and will write a generax mapping file')
parser.add_argument("-f","--fasta",
					type=str,
					help="Fasta alignment")
parser.add_argument("-s","--starting_tree",
					type=str,
					default='',
					help="Starting tree")
parser.add_argument("-m","--model",
					type=str,
					help="Model for genetree estimation.")
parser.add_argument("-mp","--map",
					type=str,
					default='',
					help="Map file for genes and species.")
parser.add_argument("-o","--output_file",
					type=str,
					default='output',
					help="name of outputfile")
args=parser.parse_args()

fasta = args.fasta
starting_tree = args.starting_tree
model = args.model
map = args.map
output = args.output_file


# Write nexus file
outFile = open(output, "w")
outFile.write('[FAMILIES]\n')
outFile.write('-family_1\n')
if starting_tree != '':
	outFile.write('starting_gene_tree = ' + starting_tree + '\n')
outFile.write('alignment = ' + fasta + '\n')
if map != '':
	outFile.write('mapping = ' + map + '\n')
outFile.write('subst_model = ' + model + '\n')


outFile.close
		