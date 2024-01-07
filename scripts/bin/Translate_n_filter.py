#!/usr/bin/env python

import argparse
import sys, os, shutil
import subprocess
from Bio import SeqIO
from Bio.Seq import Seq
import copy

# Command line options
parser = argparse.ArgumentParser(description='This script will take an input fasta and spit out a fasta of amaino acids (.faa). It assumes that all sequences in the fasta are coding sequences and that the first site starts the reading frame (i.e. eery sequence should start with ATG)')
parser.add_argument("-f","--fasta",
					type=str,
					help="Fasta to translate")
parser.add_argument("-o","--output_prefix",
					type=str,
					help="prefix to append to output files (internal_stops.fasta, incomplete_translation_frame.fasta)")
parser.add_argument("-hs","--have_stop",
					type=str,
                    default='True',
					help="Do your CDS include the final stop codon (True or False)? Default is 'True'.")
args=parser.parse_args()

output_prefix = args.output_prefix
have_stop = args.have_stop
#########################################################################################

if have_stop not in ['True', 'TRUE','true','t','T','False','FALSE','false','f','F']:
    raise Exception('True and False are the only acceptable inputs for have_stop')

#########################################################################################


def fix_incomplete_translation_frame(sequences):
    for seq in sequences:
        if len(seq.seq) % 3 != 0:
            tf1 = seq.seq.translate()
            tf2 = seq.seq[1:].translate()
            tf3 = seq.seq[2:].translate()
            stop_counts = [tf1.count('*'), tf2.count('*'),tf3.count('*')]
            if tf1.count('*') == min(stop_counts):
                translationframe = tf1
                if min(stop_counts) > 0:
                    stop_pos = tf1.index('*')
                    codon_stop = stop_pos*3
                    seq.seq = seq.seq[0:codon_stop+1]
            elif tf2.count('*') == min(stop_counts):
                translationframe = tf2
                if min(stop_counts) > 0:
                    stop_pos = tf2.index('*')
                    codon_stop = stop_pos*3
                    seq.seq = seq.seq[1:codon_stop+2]
            elif tf3.count('*') == min(stop_counts):
                translationframe = tf3
                if min(stop_counts) > 0:
                    stop_pos = tf2.index('*')
                    codon_stop = stop_pos*3
                    seq.seq = seq.seq[2:codon_stop+3]
    return(sequences)


def find_no_stop(sequences):
    no_stop_names = [seq.id for seq in sequences if '*' not in seq.seq.translate()]
    return(no_stop_names)
    

def remove_final_stop(sequences):
    for seq in sequences:
        test_seq = copy.deepcopy(seq.seq.translate())
        if len(test_seq) > 0: 
            #print(test_seq[-1])
            #print(1)
            if test_seq[-1] == "*":
                #print(2)
                seq.seq = seq.seq[:-3]
                #print(seq.seq)
        else:
            print(seq.id)
    return(sequences)


def find_internal_stop(sequences):
    ## this function assumes final stop codons have been removed
    internal_stop_names = [seq.id for seq in sequences if '*' in seq.seq.translate()]
    return(internal_stop_names)

#########################################################################################

## Parse sequences, store as list
sequences = list(SeqIO.parse(args.fasta,"fasta"))

## First exclude incomplete frames as those are problematic in one form or another
test_sequences = fix_incomplete_translation_frame(sequences)
#test_sequences = [seq for seq in sequences if seq.id not in incomplete_frames_names]

## Next, we will check for internal stop codons, and missing stop codons if appropriate
if have_stop in ['True', 'TRUE','true','t','T']:
    ## ID and remove sequences without a stop codon
    #no_stop_names = find_no_stop(test_sequences)
    #test_sequences = [seq for seq in test_sequences if seq.id not in no_stop_names]
    ## remove final stop
    test_sequences = remove_final_stop(test_sequences)

elif  have_stop in ['False','FALSE','false','f','F']:
    pass

else:
    raise Expection("Well I don\'t know why the script ran this far but True and False are the only acceptable inputs for have_stop")


## Check for internal stop codons
internal_stop_names = find_internal_stop(test_sequences)
test_sequences = [seq for seq in test_sequences if seq.id not in internal_stop_names]

test_seq_names = [seq.id for seq in test_sequences]


## Finally, filter the original sequences
good_sequences = [seq for seq in sequences if seq.id in test_seq_names]
good_output = output_prefix + '_cds_filtered.fasta'
handle = open(good_output,"w")
SeqIO.write(good_sequences, good_output, "fasta-2line")
handle.close


#incomplete_frames = [seq for seq in sequences if seq.id in incomplete_frames_names]
#incomplete_output = output_prefix + '_incomplete_translation_frame.fasta'
#handle = open(incomplete_output,"w")
#SeqIO.write(incomplete_frames, incomplete_output, "fasta-2line")
#handle.close

internal_stops = [seq for seq in sequences if seq.id in internal_stop_names]
if len(internal_stops) > 0:
    internal_stops_output = output_prefix + '_internal_stops.fasta'
    handle = open(internal_stops_output,"w")
    SeqIO.write(internal_stops, internal_stops_output, "fasta-2line")
    handle.close


#if no_stop_names:
#    no_stops = [seq for seq in sequences if seq.id in no_stop_names]
#    no_stops_output = output_prefix + '_no_stop.fasta'
#    handle = open(no_stops_output,"w")
#    SeqIO.write(no_stops, no_stops_output, "fasta-2line")


for seq in good_sequences:
        seq.seq = seq.seq.translate()

	
name = args.fasta.split('.')[0]
file = name + '.faa'
handle = open(file,"w")
SeqIO.write(good_sequences, file, "fasta-2line")
handle.close

		
## Write results to file
## Line below hashed out in V2 to write sequences as a single line


