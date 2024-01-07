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
from Bio.Seq import MutableSeq
from Bio.SeqIO import FastaIO
from dfply import *
from ete3 import Tree
from ete3 import PhyloTree

########################################
############### ARGUMENTS ##############
########################################


parser = argparse.ArgumentParser(description='')
#parser.add_argument("-og","--orthogroup",
#					type=str,
#					help="Orthogroup from an OrthoFinder orthogroups.tsv file")
parser.add_argument("-o","--orthogroup",
					type=str,
					help="Name of the orthogroup")
parser.add_argument("-a","--alignments_dir",
					type=str,
                    default = "3.3_OG_alignments",
					help="Directory containing fasta seqence alignment. Will be converted to phylip format")
parser.add_argument("-oc","--orthogroup_classifications_dir",
					type=str,
                    default="",
					help=". Default assumes it is in path")
parser.add_argument("-cd","--codeml_dir",
					type=str,
                    default="4.2_OG_codeml",
					help="name of codeml output dir")
parser.add_argument("-sd","--scripts_dir",
					type=str,
                    default="../Comparative_genomics_snakes.git/scripts/bin",
					help="Directory with scripts for the pipeline")
parser.add_argument("-cp","--codeml_path",
					type=str,
                    default="codeml",
					help="full path to codeml. Default assumes it is already in your path.")
parser.add_argument("-m","--mafft",
					type=str,
                    default="mafft",
					help="full path to mafft. Default assumes it is already in your path.")
args = parser.parse_args()

orthogroup = args.orthogroup
alignments_dir = args.alignments_dir
orthogroup_classifications_dir = args.orthogroup_classifications_dir
codeml_dir = args.codeml_dir
scripts_dir = args.scripts_dir
codeml_path = args.codeml_path
mafft = args.mafft

########################################
################# SETUP ################
########################################
nucleotides = ['A','C','G','T','-']



########################################
############## FUNCTIONS ###############
########################################

def filter_align_seqs(sequences, tree, mafft):
    current_dir = os.getcwd()
    sample_dir = codeml_dir + "/" + orthogroup
    os.chdir(sample_dir)
    focal_seqs_names = tree.get_leaf_names()
    focal_seqs = [seq for seq in sequences if seq.id in focal_seqs_names]
    for seq in focal_seqs:
        seq.seq = seq.seq.ungap()
    file = "tmp.fna"
    handle = open(file, "w")
    SeqIO.write(focal_seqs, file, "fasta")
    handle.close()
    for seq in focal_seqs:
        seq.seq=seq.seq.translate()
    file = "tmp.faa"
    handle = open(file, "w")
    SeqIO.write(focal_seqs, file, "fasta")
    handle.close()
    mafft_call = mafft + " --leavegappyregion tmp.faa > tmp_aln.faa"
    sp.call(mafft_call, shell = True)
    pal2nal_call =  "pal2nal.pl tmp_aln.faa tmp.fna -output fasta > tmp_aln.fna"
    sp.call(pal2nal_call, shell = True)
    final_seqs = list(SeqIO.parse("tmp_aln.fna","fasta"))
    for seq in final_seqs:
        seq_str = MutableSeq(seq.seq)
        for site in range(0,len(seq_str)):
            if seq_str[site] not in nucleotides:
                seq_str[site] = 'N'
        seq.seq = seq_str
    file = orthogroup + ".phy"
    out = open(file, 'w')
    out.write(' '+str(len(final_seqs))+' '+str(len(final_seqs[0].seq)))
    for seq in final_seqs:
        out.write("\n")
        out.write(seq.id+"\n")
        string=str(seq.seq)
        out.write(string)
    out.close()
    os.chdir(current_dir)


def mod_tree(tree):
    for leaf in tree.get_leaves():
        new_name = leaf.name.split('_',1)[1]
        leaf.name = new_name
    outfile = codeml_dir + "/" + orthogroup + "/" + orthogroup + ".nwk"
    tree.write(format=5,outfile=outfile)
    return(tree)


def make_codeml_ctl(orthogroup,codeml_dir,codeml_path):
    mycwd = os.getcwd()
    codeml_sample = codeml_dir + "/" + orthogroup
    os.chdir(codeml_sample)
    codeml_ctl = "python ../../" + scripts_dir + "/Format_CodeML_ctl.py -s " +  orthogroup + ".phy -t " + orthogroup + ".nwk -o " + orthogroup + "_output -r 0 -st 1 -cf 2 -nd 1 -c 0 -ad 0 -ar 0 -m 1 -ns 0 -i 0 -mg 0 -fk 0 -om 0.4 -fa 1 -a 0 -ncg 8 -gse 0 -ra 1 -cd 0 -met 0"
    sp.call(codeml_ctl, shell=True)
    codeml_call = codeml_path + " codeml.ctl > codeml.log"
    sp.call(codeml_call, shell=True)
    os.chdir(mycwd)


########################################
################# CODE #################
########################################

seq_file = alignments_dir + "/" + orthogroup + "/" + orthogroup + "_aln.fasta"
sequences = list(SeqIO.parse(seq_file,"fasta"))

tree_file = orthogroup_classifications_dir + "/" + orthogroup + "/" + orthogroup + "_codeml.nwk"
tree = PhyloTree(tree_file)

tree = mod_tree(tree)

filter_align_seqs(sequences, tree, mafft)

make_codeml_ctl(orthogroup,codeml_dir,codeml_path)



