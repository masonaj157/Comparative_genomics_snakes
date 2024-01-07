#!/usr/bin/env python

import argparse
import csv
import pandas as pd
import sys, os, shutil
import subprocess as sp
from scipy.stats import gmean
import math
import numpy as np
from skbio.stats.composition import multiplicative_replacement


# Command line options
parser = argparse.ArgumentParser(description='This script takes the output of feature and creates a csv file formatted for A. J. Mason\'s Rscripts to make venom gland transcriptome figures')
parser.add_argument("-og","--orthogroup",
					type=str,
					help="Orthogroup you want to collect read mapping info from.")
parser.add_argument("-ed","--expression_dir",
					type=str,
					default='1_read_mapping',
					help="Directory containing species-specific directories \(named after species\) containing output files from featureCounts for each species. Default is '1_read_mapping'")
parser.add_argument("-gd","--genomes_dir",
					type=str,
					default='0_Inputs/genomes',
					help="directory containing genome files for the focal species. Used to ID focal species")
parser.add_argument("-oc","--orthogroup_classifications_dir",
					type=str,
					default='3.6_OG_classifications',
					help="directory containing the paralog classifications for an orthogroup")
parser.add_argument("-ct","--counts_tpm",
					type=str,
                    default='tpm',
					help="string that identifies whether you would like your ouput as read counts or transcripts per million (tpm). Options are 'count','tpm','tpm10k','clrTPM','clrTPM10K'. TPM is default")
parser.add_argument("-o","--output_dir",
					type=str,
                    default="3.8_OG_expression",
                    help="directory of where to output the reformatted expression table")                    
args=parser.parse_args()


#Toxins = ['3FTx','BPP','CRISP','CTL','DIS','Ficolin','FusedToxin','HYAL','KAZ','KUN','LAO','LAAO','LIPA','MYO','NGF','NUC','PDE','PLA2','PLB','SNACLEC','SVMPI','SVMPII','SVMPIII','SVMP','SVSP','TruncHYAL','VEGF','Vespryn','VF','Waprin','VDP4','GluCyc']


orthogroup = args.orthogroup
expression_dir = args.expression_dir
genomes_dir = args.genomes_dir
orthogroup_classifications_dir = args.orthogroup_classifications_dir
counts_tpm = args.counts_tpm
output_dir = args.output_dir



#featurecounts_file = args.featurecounts_files
#species_list = args.species_list
#orthogroups_file = args.orthogroups_file
#output = args.output
#########################################################################################


def find_species(genomes_dir) :
	species_list = list(os.listdir(genomes_dir))
	species_list = [x.split('.')[0] for x in species_list]
	return(species_list)


def collect_feature_counts(focal_species):
    species_feature_counts = []
    for species in focal_species:
        counts_file = expression_dir + "/" + species + "/" + species + "_counts"
        counts = pd.read_csv(counts_file , skiprows=1, sep='\t')
        start_col = counts.columns.get_loc("Geneid")
        counts = counts.iloc[:,start_col:]
        species_feature_counts.append(counts)
    return(species_feature_counts)


def counts_to_TPM(species_feature_counts):
    for species in species_feature_counts:
        #print(species.columns)
        for sample in species.columns[6:]:
            Lengths = [float(x / 1000) for x in species["Length"]]
            counts = [float(x) for x in species[sample]]
            RPK = [float(count / length) for count, length in zip(counts,Lengths)]
            scaling_factor = float(sum(RPK) / 1000000)
            if scaling_factor == 0:
                scaling_factor += 1
                print('huh, you must not have any expression for one of these samples')
            TPM = [float(RPK_val / float(scaling_factor)) for RPK_val in RPK]
            species[sample] = TPM
    return(species_feature_counts)


def counts_to_clrTPM(species_feature_counts):
    for species in species_feature_counts:
        #print(species.columns)
        for sample in species.columns[6:]:
            Lengths = [float(x / 1000) for x in species["Length"]]
            counts = [float(x) for x in species[sample]]
            RPK = [float(count / length) for count, length in zip(counts,Lengths)]
            scaling_factor = float(sum(RPK) / 1000000)
            if scaling_factor == 0:
                scaling_factor += 1
                print('huh, you must not have any expression for one of these samples')
            TPM = [float(RPK_val / float(scaling_factor)) for RPK_val in RPK]
            TPM = multiplicative_replacement(np.array(TPM)).tolist()
            #print(TPM)
            geo_mean = float(gmean(TPM))
            #print(geo_mean)
            clrTPM = [float(math.log(float(TPM_val / geo_mean))) for TPM_val in TPM]
            species[sample] = clrTPM
    return(species_feature_counts)


def counts_to_TPM10K(species_feature_counts):
    for species in species_feature_counts:
        #print(species.columns)
        for sample in species.columns[6:]:
            Lengths = [float(x / 1000) for x in species["Length"]]
            counts = [float(x) for x in species[sample]]
            RPK = [float(count / length) for count, length in zip(counts,Lengths)]
            scaling_factor = float(sum(RPK) / 1000000)
            if scaling_factor == 0:
                scaling_factor += 1
                print('huh, you must not have any expression for one of these samples')
            TPM = [float(RPK_val / float(scaling_factor)) for RPK_val in RPK]
            print(len(TPM))
            TPM10K = [float((TPM_val * float(len(TPM))) / float(10000)) for TPM_val in TPM ]
            species[sample] = TPM10K
    return(species_feature_counts)


def counts_to_clrTPM10K(species_feature_counts):
    for species in species_feature_counts:
        #print(species.columns)
        for sample in species.columns[6:]:
            Lengths = [float(x / 1000) for x in species["Length"]]
            counts = [float(x) for x in species[sample]]
            RPK = [float(count / length) for count, length in zip(counts,Lengths)]
            scaling_factor = float(sum(RPK) / 1000000)
            if scaling_factor == 0:
                scaling_factor += 1
                print('huh, you must not have any expression for one of these samples')
            TPM = [float(RPK_val / float(scaling_factor)) for RPK_val in RPK]
            #print(len(TPM))
            TPM10K = [float((TPM_val * float(len(TPM))) / float(10000)) for TPM_val in TPM ]
            TPM10K = multiplicative_replacement(np.array(TPM10K)).tolist()
            geo_mean = float(gmean(TPM10K))
            clrTPM10K = [float(math.log(float(TPM10K_val / geo_mean))) for TPM10K_val in TPM10K]
            species[sample] = clrTPM10K
    return(species_feature_counts)


def format_ortholog_expression(species_feature_counts,ortho_keys):
    output_table = []
    header = ["Orthogroup"]
    for species in species_feature_counts:
        for sample in list(species.columns[6:]):
            sample_name = sample.split('/')[-1].split('_aln')[0]
            header.append(sample_name)
    output_table.append(header)
    for ortho in list(sorted(set(ortho_keys["orthogroup"]))):
        entry = [ortho]
        relevant_seqs = ortho_keys[ortho_keys["orthogroup"] == ortho]
        for species in species_feature_counts:
            #print(list(species["Geneid"]))
            #print('Smili-SVMPo-RA' in list(species["Geneid"]))
            #if not any(list(species["Geneid"]) in list(relevant_seqs["gene"])):
            zeros_entry = []
            for col in species.columns[6:]:
                #print(min(list(species[col])))
                zeros_entry.append(min(list(species[col])))
            if not any([ x in list(species["Geneid"]) for x in list(relevant_seqs["gene"])]):
                #zeros_entry = [0] * len(list(species.columns[:]))
                [entry.append(x) for x in zeros_entry]
            else:
                seqname = [ x for x in list(species["Geneid"]) if x in list(relevant_seqs["gene"]) ]
                if len(seqname) != 1:
                    print("WTF?")
                    print(seqname)
                else:
                    match = seqname[0]
                    seqcounts = species[species["Geneid"] == match].values.flatten().tolist()[6:]
                    [entry.append(x) for x in seqcounts]
        output_table.append(entry)
    return(output_table)
        

#########################################################################################

focal_species = find_species(genomes_dir)
if '' in focal_species:
    focal_species.remove('')
if '.dir' in focal_species:
    focal_species.remove('.dir')


species_feature_counts = collect_feature_counts(focal_species)


ortho_keys_file = orthogroup_classifications_dir + "/" + orthogroup + "/" + orthogroup + "_genes.csv"
ortho_keys = pd.read_csv(ortho_keys_file, sep=',')


if counts_tpm == 'tpm' or counts_tpm == 'TPM':
    species_feature_counts = counts_to_TPM(species_feature_counts)
elif counts_tpm == 'tpm10k' or counts_tpm == 'TPM10K':
    species_feature_counts = counts_to_TPM10K(species_feature_counts)
elif counts_tpm == 'clrtpm' or counts_tpm == 'clrTPM' or counts_tpm == 'CLRTPM':
    species_feature_counts = counts_to_clrTPM(species_feature_counts)
elif counts_tpm == 'clrtpm10k' or counts_tpm == 'clrTPM10K' or counts_tpm == 'CLRTPM10K' or counts_tpm == 'clrTPM10k':
    species_feature_counts = counts_to_clrTPM10K(species_feature_counts)


ortho_expression = format_ortholog_expression(species_feature_counts,ortho_keys)


mkdir_call = "mkdir " + output_dir + "/" + orthogroup
sp.call(mkdir_call, shell=True)

outfile = output_dir + "/" + orthogroup + "/" + orthogroup + "_" + counts_tpm + "_expression.csv"
with open(outfile, 'w') as csv_file:
    csv_writer = csv.writer(csv_file, delimiter = ',')
    for row in ortho_expression:
        csv_writer.writerow(row)		

csv_file.close()


