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
#parser.add_argument("-r","--read_mapping_dir",
#					type=str,
#                    default='1_read_mapping',
#					help="Directory contining gtfs that were used in read mapping. Will be used to extract gene position from transcripts")
parser.add_argument("-c","--orthogroup_classifications",
					type=str,
					default='3.6_OG_classifications',
					help="name of the output CDS file")
parser.add_argument("-go","--gene_order_dir",
					type=str,
					default='3.9_OG_geneorder',
					help="name of the output CDS file")
parser.add_argument("-d","--decostar_dir",
					type=str,
					default='3.10_OG_decostar',
					help="directory of where to write the input files for decostar")
args=parser.parse_args()

orthogroup = args.orthogroup
genomes_dir = args.genomes_dir
gene_order_dir = args.gene_order_dir
decostar_dir = args.decostar_dir
#read_mapping_dir = args.read_mapping_dir
orthogroup_classifications = args.orthogroup_classifications



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



def make_adjacencies(gene_order,focal_species):
    output_table = []
    for species in focal_species:
        subset = gene_order[gene_order["species"] == species]
        subset = subset.sort_values(by=['seq','start'])
        print(subset["original_id"])
        print(subset["start"])
        for i in range(1,len(subset)):
            entry1 = species + '_' + subset["original_id"].iloc[i-1]
            entry2 = species + '_' + subset["original_id"].iloc[i]
            table_row =[entry1,entry2]
            output_table.append(table_row)
    return(output_table)

########################################
################# CODE #################
########################################

focal_species = find_species(genomes_dir)
if '' in focal_species:
    focal_species.remove('')
if '.dir' in focal_species:
    focal_species.remove('.dir')


gene_order_file = gene_order_dir + "/" + orthogroup + "/" + orthogroup + "_positions.csv"
gene_order = pd.read_csv(gene_order_file, ",")

output_table = make_adjacencies(gene_order, focal_species)

make_output_folder = "mkdir " + decostar_dir + '/' + orthogroup
subprocess.call(make_output_folder, shell=True)

## Write the decostar adjacencies file
adjacencies_table_outfile = decostar_dir + '/' + orthogroup + "/" + orthogroup + '_adjacencies.txt'
with open(adjacencies_table_outfile, 'w') as csv_file:
    csv_writer = csv.writer(csv_file, delimiter = '\t')
    for row in output_table:
        csv_writer.writerow(row)		
csv_file.close()



make_tree_folder = "mkdir " + decostar_dir + '/' + orthogroup + "/" + orthogroup + "_trees"
subprocess.call(make_tree_folder, shell=True)

copy_gene_tree = "cp " + orthogroup_classifications + "/" + orthogroup + "/" + orthogroup + "_decostar.nwk " + decostar_dir + '/' + orthogroup + "/" + orthogroup + "_trees/"
copy_species_tree = "cp " + orthogroup_classifications + "/" + orthogroup + "/" + orthogroup + "_pruned_species.nwk " + decostar_dir + '/' + orthogroup + "/"
subprocess.call(copy_gene_tree, shell=True)
subprocess.call(copy_species_tree, shell=True)



## write the decostar distrib_gene_trees.txt file
distrib_gene_trees_file = decostar_dir + '/' + orthogroup + "/" + orthogroup + "_distrib_gene_trees.txt"
outFile = open(distrib_gene_trees_file, "w")
outFile.write(orthogroup + "_trees/" + orthogroup + "_decostar.nwk")
outFile.close

param_file = decostar_dir + '/' + orthogroup + "/" + orthogroup + ".param.txt"
outFile = open(param_file, "w")
outFile.write('species.file=' + orthogroup + "_pruned_species.nwk" + '\n')
outFile.write('gene.distribution.file=' + orthogroup + "_distrib_gene_trees.txt" + '\n')
outFile.write('adjacencies.file=' + orthogroup + '_adjacencies.txt' + '\n')
outFile.write('with.transfer=false' + '\n')
outFile.write('\n')
outFile.write('verbose=1' + '\n')
#outFile.write('rooted=true' + '\n')
outFile.write('superverbose=0' + '\n')
outFile.write('write.newick=0' + '\n')
outFile.write('\n')
outFile.write('dated.species.tree=0' + '\n')
outFile.write('\n')
outFile.write('char.sep=_' + '\n')
outFile.write('\n')
#outFile.write('try.all.amalgamation=0' + '\n')
outFile.write('\n')
#outFile.write('ARt-DeCo=1' + '\n')
outFile.write('write.newick=true' + '\n')
outFile.write('write.genes=true' + '\n')
outFile.write('write.adjaceny.trees=true' + '\n')
outFile.write('\n')
outFile.write('output.dir=output')
outFile.close
