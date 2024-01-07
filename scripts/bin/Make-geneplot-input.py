#!/usr/bin/env python


import argparse
import csv
import pandas as pd
import os
import subprocess
from Bio import SeqIO
from Bio import Seq
from Bio.SeqRecord import SeqRecord

parser = argparse.ArgumentParser(description='A script to extract gene or CDS sequences from a genome given a genome fasta and gff. Sadly quite slow, but I\'m no computational biologist. ')
parser.add_argument("-o","--orthogroup",
					type=str,
					default='',
					help="orthogroup")
parser.add_argument("-g","--genomes_dir",
					type=str,
                    default='0_Inputs/genomes',
					help='Directory containing input genomes. Only used to get species names')
parser.add_argument("-r","--read_mapping_dir",
					type=str,
                    default='1_read_mapping',
					help="Directory contining gtfs that were used in read mapping. Will be used to extract gene position from transcripts")
parser.add_argument("-c","--orthogroup_classifications",
					type=str,
					default='3.6_OG_classifications',
					help="name of the output CDS file")
parser.add_argument("-out","--out_dir",
					type=str,
					default='3.9_OG_geneorder',
					help="name of the output CDS file")
args=parser.parse_args()

orthogroup = args.orthogroup
genomes_dir = args.genomes_dir
read_mapping_dir = args.read_mapping_dir
orthogroup_classifications = args.orthogroup_classifications
out_dir = args.out_dir

########################################
############## FUNCTIONS ###############
########################################

def find_species(genomes_dir) :
	species_list = list(os.listdir(genomes_dir))
	species_list = [x.split('.')[0] for x in species_list]
	return(species_list)


def collect_feature_counts(focal_species):
    species_feature_counts = []
    for species in focal_species:
        counts_file = read_mapping_dir + "/" + species + "/" + species + "_counts"
        counts = pd.read_csv(counts_file , skiprows=1, sep='\t')
        species_feature_counts.append(counts)
    return(species_feature_counts)



def make_table(focal_species,counts_list,gene_info):
    output_table = []
    for species,species_counts in zip(focal_species,counts_list):
        for gene in list(gene_info["gene"]):
            relevant_entry = species_counts[species_counts["Geneid"] == gene]
            if len(relevant_entry) == 1:
                gene_entry = list(gene_info[gene_info["gene"] == gene]["orthogroup"])[0]
                start = int(relevant_entry["Start"].iloc[0].split(';')[0])
                end = int(relevant_entry["End"].iloc[0].split(';')[-1])
                strand_sign = relevant_entry["Strand"].iloc[0].split(';')[0]
                if strand_sign == '+':
                    strand = "forward"
                    orientation = 1
                elif strand_sign == "-":
                    strand = "reverse"
                    orientation = -1
                seq = relevant_entry["Chr"].iloc[0].split(';')[0]
                molecule = species + '-' + seq
                original_id = gene
                entry = [molecule,species,gene_entry,start,end,strand,orientation,seq,original_id]
                output_table.append(entry)
    return(output_table)

########################################
################# CODE #################
########################################

focal_species = find_species(genomes_dir)
if '' in focal_species:
    focal_species.remove('')
if '.dir' in focal_species:
    focal_species.remove('.dir')


counts_list = collect_feature_counts(focal_species)

gene_info_file = orthogroup_classifications + "/" + orthogroup + "/" + orthogroup + "_genes.csv"
gene_info = pd.read_csv(gene_info_file, ",")

output_table = make_table(focal_species,counts_list,gene_info)


make_output_folder = "mkdir " + out_dir + '/' + orthogroup
subprocess.call(make_output_folder, shell=True)

classes_table_outfile = out_dir + '/' + orthogroup + "/" + orthogroup + '_positions.csv'
with open(classes_table_outfile, 'w') as csv_file:
    csv_writer = csv.writer(csv_file, delimiter = ',')
    csv_writer.writerow(['molecule','species','gene','start','end','strand','orientation','seq','original_id'])
    for row in output_table:
        csv_writer.writerow(row)		
csv_file.close()





