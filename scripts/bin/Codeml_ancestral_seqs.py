#!/usr/bin/env python


import argparse
import csv
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
#from Bio.Alphabet import IUPAC

parser = argparse.ArgumentParser(description='A script to extract gene or CDS sequences from a genome given a genome fasta and gff. Sadly quite slow, but I\'m no computational biologist. ')
parser.add_argument("-r","--rst_file",
					type=str,
					default='',
					help="The output rst file from a codeml run")
parser.add_argument("-o","--output",
					type=str,
					default='output.fasta',
					help="Name of the output fasta (e.g. output.fasta)")
args=parser.parse_args()

blast_file = args.blast_file
output = args.output


################################################

def rst_parse(rst_file) :
    rst_list = []
    rst_res = open(rst_file,'r')
    Lines = rst_res.readlines()
    first = True
    for line in Lines:
        anc_seqs = []
        if "tree with node labels for Rod Page's TreeView" in line:
            read_tree = True
        if read_tree = True:
            tree_string = line.split('\n')[0]
            tree_string = [character for character in tree_string if character != ' ']
            tree_string = "".join(tree_string)
            tree = PhyloTree(tree_string)
            read_tree = False
        if "List of extant and reconstructed sequences" in line:
            found_seqs = True
            first_node = False
        if found_seqs == True and line[0:3] == "node":
            first_node = True
            name = line.split(' ')[1].split('#')[1]
            seq_string = line[18:].split("\n")[0]
            seq = [character for character in seq_string if character != ' ']
            seq = "".join(seq)
            new_seq = SeqRecord(Seq(seq),id=name)
            anc_seqs.append(new_seq)
        elif found_seq == True and line == "\n":
            break
    return(tree, anc_seqs)


################################################

sequences = rst_parse(rst_file)

handle = open(output, "w")
SeqIO.write(sequences, output, "fasta-2line")
handle.close()