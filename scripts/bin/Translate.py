#!/usr/bin/env python

import argparse
import sys, os, shutil
import subprocess
from Bio import SeqIO
from Bio.Seq import Seq

# Command line options
parser = argparse.ArgumentParser(description='This script will take an input fasta and spit out a fasta of amaino acids (.faa). It assumes that all sequences in the fasta are coding sequences and that the first site starts the reading frame (i.e. eery sequence should start with ATG)')
parser.add_argument("-f","--fasta",
					type=str,
					help="Fasta to translate")

args=parser.parse_args()

## Parse sequences, store as list
#sequences = list(SeqIO.parse("EOG0907004F.fasta","fasta"))
sequences = list(SeqIO.parse(args.fasta,"fasta"))

for seq in sequences:
        seq.seq = seq.seq.translate()

## check for internal stop
Good = [ seq for seq in sequences if '*' not in seq.seq[0:-1] ]
Bad = [ seq for seq in sequences if '*' in seq.seq[0:-1]]



## Remove final stop from Good
#for seq in Good:
#	seq.seq = seq.seq[0:-1]
	
	
name = args.fasta.split('.')[0]
file = name + '.faa'
		
## Write results to file
## Line below hashed out in V2 to write sequences as a single line

handle = open(file, "w")
SeqIO.write(Good, file, "fasta-2line")
handle.close()

if len(Bad) != 0:
	bad_file = name + '_internal_stop.faa'
	handle = open(bad_file, "w")
	SeqIO.write(Bad, bad_file, "fasta-2line")
	handle.close()


#if len(No_stop) != 0:
#	no_stop_file = name + '_no_stop.faa'
#	handle = open(no_stop_file, "w")
#	SeqIO.write(No_stop, no_stop_file, "fasta-2line")
#	handle.close()