#!/usr/bin/env python

# Additional software necessary to run this:
# (1) bwa 
# (2) samtools
# (3) bedtools
# (4) biopython
# (5) pandas
# (6) dfply

import argparse
from asyncio import subprocess
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
parser.add_argument("-fc","--feature_counts",
					type=str,
					help="Output of feature counts")
parser.add_argument("-cds","--cds",
					type=str,
					help="fasta of cds sequences")
parser.add_argument("-g","--gff",
					type=str,
					help="gff for the species")
args = parser.parse_args()

########################################
################# SETUP ################
########################################

feature_counts = args.feature_counts
cds = args.cds
gff = args.gff



print("""
Checking that featureCounts and extracted CDS names match
""")

#print("\n"+dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" ::: Starting...")
#print("\tInput --> "+ input_name)
#print("\tReads --> "+ os.path.basename(reads))
#print("\tRemoving transcripts with...")
#print("\t\t< " + str(cov) + "x coverage")
#print("\t\tfor > " + str(prop) + "% of the total transcript length")
#print("\tThreads --> " + str(threads))

########################################
############## FUNCTIONS ###############
########################################

class GFF:
    def __init__(self, gff_entry):
        self.seqid = gff_entry[0]
        self.source = gff_entry[1]
        self.type = gff_entry[2]
        self.start = int(gff_entry[3])
        self.end = int(gff_entry[4])
        self.score = gff_entry[5]
        self.strand = gff_entry[6]
        self.phase = gff_entry[7]
        if len(gff_entry) > 8:
            self.attributes = parse_GFF_attributes(gff_entry[8])
        else:
            self.attributes = parse_GFF_attributes("ID=unamed")
            print(gff_entry)
            print(gff_entry[0], gff_entry[2],gff_entry[3],gff_entry[4])
        if self.attributes.ID == []:
            self.attributes.ID = ['unnamed']


class gff_attributes:
    def __init__(self, gff_attributes):
        self.ID = gff_attributes[0]
        self.Name = gff_attributes[1]
        self.Alias = gff_attributes[2]
        self.Note = gff_attributes[3]
        self.Parent = gff_attributes[4]        



def GFF_parse(gff) :
	gff_list = []
	with open(gff) as OF:
		reader = csv.reader(OF, delimiter='\t')
		for row in reader :
			if (row[0][0] != '#') and (len(row) > 1) :
				gff_entry = GFF(row)
				gff_list.append(gff_entry)
	return(gff_list)


def parse_GFF_attributes(gff_attribute) :
    split_attributes_1 = list(gff_attribute.split(';'))
    split_attributes = []
    for attribute in split_attributes_1:
        split_attributes.append(list(attribute.split('=')))
    ID = [x[1] for x in split_attributes if x[0] == 'ID']
    Name = [x[1] for x in split_attributes if x[0] == 'Name']
    Alias = [x[1] for x in split_attributes if x[0] == 'Alias']
    Note = [x[1] for x in split_attributes if x[0] == 'Note']
    Parent = [x[1] for x in split_attributes if x[0] == 'Parent']
    attributes_list = gff_attributes([ID, Name, Alias, Note, Parent])
    return(attributes_list)



def pull_feature_names(counts) :
    feature_names = list(counts["Geneid"])
    return(feature_names)


def test_and_id_mismatches(missing_cds,  missing_features):
    seq_names = [seq.id for seq in sequences]
    #missing_cds = [x for x in seq_names if x not in feature_names]
    if len(missing_cds) > 0 :
        print(str(len(missing_cds)) + " cds without a match")
        print(missing_cds)
    #missing_features = [x for x in feature_names if x not in seq_names]
    if len(missing_features) > 0:
        print(str(len(missing_features)) + " features without a match")
        print(missing_features)
    if len(missing_cds) > 0 or len(missing_features) > 0:
        return(True)
    else:
        return(False)



def fix_bad_parents(missing_features, counts, gff):
    print("Reading gff")
    gff_list = GFF_parse(gff)
    cds_gff_entries = [entry for entry in gff_list if entry.type == 'CDS']
    for feature in missing_features:
        relevant_gff_entry = [entry for entry in cds_gff_entries if len(entry.attributes.Parent) != 0 and entry.attributes.Parent[0] == feature]
        if len(relevant_gff_entry) != 0:
            ID = relevant_gff_entry[0].attributes.ID
            counts["Geneid"] = counts['Geneid'].replace([feature],ID)
    return(counts)


def remove_cds_label_seqs(sequences):
    for seq in sequences:
        if ":cds" in seq.id:
            seq.id = seq.id.split(':cds')[0]
    return(sequences)


def remove_cds_label_features(counts):
    for entry in list(counts["Geneid"]):
        if ":cds" in entry:
            new_id = entry.split(':cds')[0]
            counts["Geneid"] = counts['Geneid'].replace([entry],new_id)
    return(counts)


def remove_unnamed_seq(sequences):
    sequences = [seq for seq in sequences if "unnamed" not in seq.id]
    return(sequences)


########################################
################# CODE #################
########################################

## Create CDS files from genomes and annotations
sequences = list(SeqIO.parse(cds,"fasta"))

counts = pd.read_csv(feature_counts, skiprows=1, sep='\t')

feature_names = pull_feature_names(counts)

seq_names = [seq.id for seq in sequences]

missing_features = [x for x in feature_names if x not in seq_names]

missing_cds = [x for x in seq_names if x not in feature_names]

if test_and_id_mismatches(missing_cds, missing_features):
    counts = fix_bad_parents(missing_features, counts, gff)
    #missing_cds = id_mismatches(sequences, feature_names)
    counts = remove_cds_label_features(counts)
    sequences = remove_cds_label_seqs(sequences)
    sequences = remove_unnamed_seq(sequences)



    feature_names = pull_feature_names(counts)

    seq_names = [seq.id for seq in sequences]

    missing_features = [x for x in feature_names if x not in seq_names]

    missing_cds = [x for x in seq_names if x not in feature_names]

    if test_and_id_mismatches(missing_cds, missing_features):
        print("well, darn it.")


    remove_seq_file = "rm " + cds
    sp.call(remove_seq_file, shell=True)


    handle = open(cds, "w")
    SeqIO.write(sequences, cds, "fasta")
    handle.close()


    tmp1 = feature_counts + '_tmp1.tsv'
    counts.to_csv(tmp1, sep="\t")

    tmp2 = feature_counts + '_tmp2.tsv'
    grab_header = "head -n 1 " + feature_counts + " > " + tmp2
    sp.call(grab_header, shell=True)
    remove_fc_file = "rm " + feature_counts
    sp.call(remove_fc_file, shell=True)
    add_count_info = "cat " + tmp1 + " >> " + tmp2
    sp.call(add_count_info, shell=True)
    rename_tmp2 = "mv " + tmp2 + " " + feature_counts
    sp.call(rename_tmp2, shell=True)
    remove_tmp1 = "rm " + tmp1
    sp.call(remove_tmp1, shell=True)


else:
    print("Everything matches!")