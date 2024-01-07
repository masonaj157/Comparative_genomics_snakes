#!/usr/bin/env python

# Additional software necessary to run this:

# (1) biopython
# (2) pandas
# (3) dfply
# (4) mafft
# (5) iqtree
# (6) generax v.XXX

import argparse
from ctypes import alignment
import sys, os, shutil
import subprocess as sp
import datetime as dt
import csv
#import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqIO import FastaIO
from dfply import *

########################################
############### ARGUMENTS ##############
########################################

parser = argparse.ArgumentParser(description='')
#parser.add_argument("-og","--orthogroup",
#					type=str,
#					help="Orthogroup from an OrthoFinder orthogroups.tsv file")
parser.add_argument("-og","--orthogroup",
					type=str,
					help="Orthogroup to process.")
parser.add_argument("-m","--species_map",
					type=str,
					help="Comma separated file file where each lines lists a species and its species assignment with headers \'Sample,Species\'")
parser.add_argument("-s","--scripts",
					type=str,
                    default="../Comparative_genomics_snakes.git/scripts/bin",
					help="Directory containing associated scripts")
parser.add_argument("-st","--species_tree",
					type=str,
                    default="0_Inputs/species_tree.nwk",
					help="Species tree file for orthofinder. Newick format")
parser.add_argument("-gd","--genomes_dir",
					type=str,
                    default="0_Inputs/genomes",
					help="Directory containing the genomes for focal taxa.")
parser.add_argument("-otn","--other_taxa_fna",
					type=str,
                    default="0_Inputs/other_taxa_fna",
					help="Directory containing the CDS nucleotide sequences for other taxa")
parser.add_argument("-ota","--other_taxa_faa",
					type=str,
                    default="0_Inputs/other_taxa_faa",
					help="Directory containing the AA sequences for other taxa")
parser.add_argument("-ed","--expression_dir",
					type=str,
                    default="1_read_mapping",
					help="Directory containing readmapping and featureCounts output for each species. Default is '1_read_mapping'")
parser.add_argument("-fna","--fna_dir",
					type=str,
                    default="2.1_CDS_fna",
					help="Directory containing the CDS nucleotide sequences for focal taxa")
parser.add_argument("-faa","--faa_dir",
					type=str,
                    default="2.2_CDS_faa",
					help="Directory containing the AA sequences for focal taxa")
parser.add_argument("-md","--map_dir",
					type=str,
                    default="3.1_Species_maps",
					help="Directory containing sequence")
parser.add_argument("-sd","--sequence_dir",
					type=str,
                    default="3.2_OG_sequences",
					help="Directory of where to write a fasta of collected orthogroup sequences")
parser.add_argument("-ad","--alignment_dir",
					type=str,
                    default="3.3_OG_alignments",
					help="Directory of where to write an aligned fasta of sequences")
parser.add_argument("-td","--tree_dir",
					type=str,
                    default="3.4_OG_trees",
					help="Directory of where to run iqtree")
parser.add_argument("-rd","--reconciliations_dir",
					type=str,
                    default="3.5_OG_reconciliations",
					help="Directory of where to run generax to reconcile the genetree")
parser.add_argument("-pcd","--paralog_classifications_dir",
					type=str,
                    default="3.6_OG_classifications",
					help="Directory of where to write the output of orthocaller which generates stats input files")
parser.add_argument("-akd","--AA_kmers_dir",
					type=str,
                    default="3.7_OG_AA_kmers",
					help="Directory of where to write the output of orthocaller which generates stats input files")
parser.add_argument("-oe","--orthogroup_expression_dir",
					type=str,
                    default="3.8_OG_expression",
					help="Directory of where to write the output of expression orthogroup formatter. Will produce an orthogroup specific directory and csv file ascribing expression to genes within an orthogroup based on classifications. Default is '3.8_OG_expression'")
parser.add_argument("-gp","--gene_positions_dir",
					type=str,
                    default="3.9_OG_geneorder",
					help="Directory of where to write the table used to visualize gene positions'")
parser.add_argument("-dod","--dissimilarity_output_dir",
					type=str,
                    default="4.1_OG_Routputs",
					help="Directory where to write dissimilarity output directory and files. Default is 4.1_OG_Routputs")
parser.add_argument("-cm","--codeml_dir",
					type=str,
                    default="4.2_OG_codeml",
					help="Directory where to write codeml output directory and files. Default is 4.2_OG_codeml")
parser.add_argument("-k","--kmer_size",
					type=int,
                    default=20,
					help="AA kmer size for counting")
#parser.add_argument("-ct","--counts_tpm",
#					type=str,
#                    default="tpm",
#					help="string that identifies whether you would like your ouput as read counts or transcripts per million (tpm). Options are 'count' or 'tpm'. TPM is default")
parser.add_argument("-mft","--mafft",
					type=str,
                    default="mafft",
					help="Path to iqtree. Default assumes it is in path")
parser.add_argument("-iq","--iqtree",
					type=str,
                    default="iqtree",
					help="Path to iqtree. Default assumes it is in path")
parser.add_argument("-cp","--codeml_path",
					type=str,
                    default="codeml",
					help="full path to codeml. Default assumes it is already in your path.")
args = parser.parse_args()


########################################
################# SETUP ################
########################################

orthogroup = args.orthogroup
genomes_dir = args.genomes_dir
species_map = args.species_map
map_dir = args.map_dir
fna_dir = args.fna_dir
faa_dir = args.faa_dir
otn_dir = args.other_taxa_fna
ota_dir = args.other_taxa_faa
sequence_dir = args.sequence_dir
alignment_dir = args.alignment_dir
tree_dir = args.tree_dir
reconciliations_dir= args.reconciliations_dir
paralog_classifications_dir = args.paralog_classifications_dir
AA_kmers_dir = args.AA_kmers_dir
dissimilarity_output_dir = args.dissimilarity_output_dir
codeml_dir = args.codeml_dir
scripts = args.scripts
species_tree = args.species_tree
expression_dir = args.expression_dir
orthogroup_expression_dir = args.orthogroup_expression_dir
gene_positions_dir = args.gene_positions_dir
kmer_size = args.kmer_size
##counts_tpm = args.counts_tpm
mafft = args.mafft
iqtree = args.iqtree
codeml_path = args.codeml_path
#out = args.output_directory


print("""
Processing output of OrthoFinder. Building alignments, trees, and reconciling.
""")

print("\n"+dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" ::: Starting...")
print("\tOrthogroup --> "+ orthogroup)
#print("\tOutput directory --> "+ out)

########################################
############## FUNCTIONS ###############
########################################

def find_species(fna_dir) :
	species_list = list(os.listdir(fna_dir))
	species_list = [x.split('.fna')[0] for x in species_list]
	return(species_list)


def pull_nuc_seqs(map,fna_dir, otn_dir) :
    seqs = []
    focal_species = find_species(fna_dir)
    for species in list(set(list(map["Species"]))):
        red_map = map[map["Species"] == species]
        if species in focal_species:
            file = fna_dir + "/" + species + ".fna"
        else:
            file = otn_dir + "/" + species + ".fna"
        fna_seqs = list(SeqIO.parse(file, "fasta"))
        for seq in fna_seqs:
            if seq.id in list(red_map["Sequence"]):
                seqs.append(seq)
    for seq in seqs:
        seq.description = ''
    return(seqs)
        

def run_mafft(mafft,orthogroup, sequence_dir,alignment_dir) :
    current_dir = os.getcwd()
    os.chdir(sequence_dir)
    translate_call = "python ../" + scripts + "/Translate.py -f " + orthogroup + ".fasta"
    sp.call(translate_call, shell = True)
    os.chdir(current_dir)
    mkdir_call = "mkdir " + alignment_dir + "/" + orthogroup
    sp.call(mkdir_call, shell=True)
    mafft_call = mafft + " --leavegappyregion " + sequence_dir + "/" + orthogroup + ".faa > " + alignment_dir + "/" + orthogroup + "/" + orthogroup + "_aln.faa"
    print(mafft_call)
    sp.call(mafft_call, shell = True)
    change_dir = alignment_dir + "/" + orthogroup
    os.chdir(change_dir)
    sp.call("pwd", shell = True)
    sp.call("ls", shell = True)
    pal2nal_call =  "pal2nal.pl " + orthogroup + "_aln.faa ../../" + sequence_dir + "/" + orthogroup + ".fasta -output fasta > " + orthogroup + "_aln.fasta"
    sp.call(pal2nal_call, shell = True)
    os.chdir(current_dir)


def run_iqtree(orthogroup, alignment_dir, tree_dir, iqtree) :
    current_dir = os.getcwd()
    mk_tree_dir = "mkdir " + tree_dir + "/" + orthogroup
    sp.call(mk_tree_dir, shell=True)
    copy_alignment = "cp " + alignment_dir + "/" + orthogroup + "/" + orthogroup + "_aln.fasta " + tree_dir + "/" + orthogroup + "/"
    sp.call(copy_alignment,shell=True)
    working_dir = tree_dir + "/" + orthogroup
    os.chdir(working_dir)
    if len(sequences) >= 4:
        iqtree_call = iqtree + " -s " + orthogroup + "_aln.fasta -m MFP -bb 1000 -nt AUTO -ntmax 8"
    else:
        iqtree_call = iqtree + " -s " + orthogroup + "_aln.fasta -m MFP -nt AUTO -ntmax 8"
    print(iqtree_call)
    sp.call(iqtree_call, shell = True)
    remove_copied_alignment = "rm " + orthogroup + "_aln.fasta"
    sp.call(remove_copied_alignment,shell=True)
    os.chdir(current_dir)


def write_generax_map(orthogroup, reconciliations_dir, map) :
    generax_map = []
    for species in list(set(list(map["Species"]))):
        red_map = map[map["Species"] == species]
        species_line = species + ":"
        for og in list(red_map["Sequence"]):
            if og == list(red_map["Sequence"])[0]:
                species_line = species_line + og
            else:
                species_line = species_line + ";" + og
        generax_map.append(species_line)
    generax_map_file = reconciliations_dir + "/" + orthogroup + "/" + orthogroup + "_generax_map.txt"
    #print(generax_map_file)
    with open(generax_map_file,'w') as outFile:
        for entry in generax_map:
                outFile.write(entry + '\n')
    outFile.close


def run_generax(orthogroup, reconciliations_dir, map, species_tree):
    current_dir = os.getcwd()
    mk_ortho_dir = "mkdir " + reconciliations_dir + "/" + orthogroup
    #print(mk_ortho_dir)
    sp.call(mk_ortho_dir, shell=True)
    #sp.call("conda activate ete3", shell=True)
    write_generax_map(orthogroup, reconciliations_dir, map)
    format_families_call = "python " + scripts + "/Format_generax_families.py -f " + current_dir + "/" + alignment_dir + "/" + orthogroup + "/" + orthogroup + "_aln.fasta -s " + current_dir + "/" + tree_dir + "/" + orthogroup + "/" + orthogroup + "_aln.fasta.treefile -m GTR+G -mp " + current_dir + "/" + reconciliations_dir + "/" + orthogroup + "/" + orthogroup + "_generax_map.txt -o " + reconciliations_dir + "/" + orthogroup + "/" + orthogroup + "_families.txt"
    sp.call(format_families_call, shell=True)
    Prune_species_tree_call = "python " + scripts + "/Prune_species_tree_by_map.py -st " + species_tree + " -m " + reconciliations_dir + "/" + orthogroup + "/" + orthogroup + "_generax_map.txt -o " + reconciliations_dir + "/" + orthogroup + "/" + orthogroup + "_species_tree.nwk"
    print(Prune_species_tree_call)
    sp.call(Prune_species_tree_call, shell=True)
    #sp.call("conda deactivate")
    working_dir = reconciliations_dir + "/" + orthogroup
    os.chdir(working_dir)
    generax_call = "generax --unrooted-gene-tree --per-species-rates -f " + orthogroup + "_families.txt -r UndatedDL --support-threshold 80 --max-spr-radius 10 -s "+ orthogroup + "_species_tree.nwk -p " + orthogroup + "_output"
    print(generax_call)
    sp.call(generax_call, shell=True)
    os.chdir(current_dir)


def run_ortho_caller(scripts, orthogroup,species_tree,reconciliations_dir, map_dir, genomes_dir, paralog_classifications_dir):
    ortho_caller_call = "python " + scripts + "/ortho_caller.py -og " + orthogroup + " -st " + species_tree + " -rd " + reconciliations_dir + " -md " + map_dir + " -gd " + genomes_dir + " -od " + paralog_classifications_dir
    sp.call(ortho_caller_call, shell=True)


def run_filter_and_collect_kmers(scripts, orthogroup, genomes_dir, aminoacid_dir, map_dir, kmer_size, output_dir):
    kmers_call = "python " + scripts + "/Filter_and_collect_AA_kmers.py -og " + orthogroup + " -gd " + genomes_dir + " -ad " + aminoacid_dir + " -m " + map_dir + " -k " + str(kmer_size) + " -o " + output_dir
    sp.call(kmers_call, shell=True)


def run_ortho_expression_formatter(scripts, orthogroup, expression_dir, genomes_dir, orthogroup_classifications_dir, outputdir):
    ortho_expression_format_call = "python " + scripts + "/Ortho_Expression_formatter_v0.1.py -og " + orthogroup + " -ed " + expression_dir + " -gd " + genomes_dir + " -oc " + orthogroup_classifications_dir + " -ct count -o " + outputdir
    sp.call(ortho_expression_format_call, shell=True)
    ortho_expression_format_call = "python " + scripts + "/Ortho_Expression_formatter_v0.1.py -og " + orthogroup + " -ed " + expression_dir + " -gd " + genomes_dir + " -oc " + orthogroup_classifications_dir + " -ct tpm -o " + outputdir
    sp.call(ortho_expression_format_call, shell=True)
    ortho_expression_format_call = "python " + scripts + "/Ortho_Expression_formatter_v0.1.py -og " + orthogroup + " -ed " + expression_dir + " -gd " + genomes_dir + " -oc " + orthogroup_classifications_dir + " -ct tpm10k -o " + outputdir
    sp.call(ortho_expression_format_call, shell=True)
    ortho_expression_format_call = "python " + scripts + "/Ortho_Expression_formatter_v0.1.py -og " + orthogroup + " -ed " + expression_dir + " -gd " + genomes_dir + " -oc " + orthogroup_classifications_dir + " -ct clrtpm -o " + outputdir
    sp.call(ortho_expression_format_call, shell=True)
    ortho_expression_format_call = "python " + scripts + "/Ortho_Expression_formatter_v0.1.py -og " + orthogroup + " -ed " + expression_dir + " -gd " + genomes_dir + " -oc " + orthogroup_classifications_dir + " -ct clrtpm10k -o " + outputdir
    sp.call(ortho_expression_format_call, shell=True)


def make_figure_inputs(scripts, orthogroup, genomes_dir, expression_dir, orthogroup_classifications_dir , gene_positions_dir):
    make_geneplot_input_call = "python " + scripts + "/Make-geneplot-input.py -o " + orthogroup + " -g " + genomes_dir + " -r " + expression_dir + " -c " + orthogroup_classifications_dir + " -out " + gene_positions_dir
    sp.call(make_geneplot_input_call,shell=True)


def run_dissimilarity_calculations(scripts,orthogroup, species_map, aminoacid_dir, orthogroups_info_dir, orthogroup_expression_dir, dissimilarity_output_dir):
    mkdir_call = "mkdir " + dissimilarity_output_dir + "/" + orthogroup
    sp.call(mkdir_call, shell=True)
    R_dissimilarity_call = "Rscript " + scripts + "/calculate_dissimilarity_v2.R -o " + orthogroup + " -a " + aminoacid_dir + " --orthogroups_info_dir " + orthogroups_info_dir + " -e " + orthogroup_expression_dir + " --output_dir " + dissimilarity_output_dir
    print(R_dissimilarity_call)
    sp.call(R_dissimilarity_call, shell=True)


def run_codeml(scripts,orthogroup,alignment_dir,paralog_classifications_dir, codeml_dir, codeml_path):
    mkdir_call = "mkdir " + codeml_dir + "/" + orthogroup
    sp.call(mkdir_call, shell=True)
    codeml_call = "python " + scripts + "/Prepare_CodeML.py -o " + orthogroup + " -a " + alignment_dir + " -oc " + paralog_classifications_dir + " -cd " + codeml_dir + " -sd " + scripts + " -cp " + codeml_path + " -m " + mafft
    print(codeml_call)
    sp.call(codeml_call, shell=True) 


########################################
################# CODE #################
########################################

map_file = map_dir + "/" + orthogroup + '.map'
map = pd.read_csv(map_file,",")


sequences = pull_nuc_seqs(map,fna_dir, otn_dir)
print(len(sequences))

seq_file = sequence_dir + "/" + orthogroup + ".fasta"
handle = open(seq_file, "w")
SeqIO.write(sequences, seq_file, "fasta-2line")
handle.close()

run_mafft(mafft,orthogroup, sequence_dir, alignment_dir)

run_iqtree(orthogroup, alignment_dir, tree_dir, iqtree)

run_generax(orthogroup, reconciliations_dir, map, species_tree)

run_ortho_caller(scripts, orthogroup,species_tree,reconciliations_dir, map_dir, genomes_dir, paralog_classifications_dir)

run_filter_and_collect_kmers(scripts, orthogroup, genomes_dir, faa_dir, map_dir, kmer_size, AA_kmers_dir)

run_ortho_expression_formatter(scripts, orthogroup, expression_dir, genomes_dir, paralog_classifications_dir, orthogroup_expression_dir)

make_figure_inputs(scripts, orthogroup, genomes_dir, expression_dir, paralog_classifications_dir, gene_positions_dir)

run_dissimilarity_calculations(scripts,orthogroup, species_map, AA_kmers_dir, paralog_classifications_dir, orthogroup_expression_dir, dissimilarity_output_dir)

run_codeml(scripts,orthogroup,alignment_dir,paralog_classifications_dir, codeml_dir, codeml_path)