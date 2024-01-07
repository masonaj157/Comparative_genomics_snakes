#!/usr/bin/env python


import argparse
from audioop import reverse
from copy import deepcopy
import csv
from email.utils import encode_rfc2231
from enum import unique
from random import seed
from typing import final
import pandas as pd
pd.options.mode.chained_assignment = None
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
parser.add_argument("-d","--decostar_dir",
					type=str,
					default='3.10_OG_decostar',
					help="directory of where to write the input files for decostar")
parser.add_argument("-c","--orthogroup_classifications",
					type=str,
					default='3.6_OG_classifications',
					help="name of the output CDS file")
parser.add_argument("-go","--geneorder_dir",
					type=str,
					default='3.9_OG_geneorder',
					help="name of the output CDS file")
args=parser.parse_args()

orthogroup = args.orthogroup
decostar_dir = args.decostar_dir
orthogroup_classifications = args.orthogroup_classifications
geneorder_dir = args.geneorder

########################################
############## FUNCTIONS ###############
########################################

def read_decostar_species(species_file):
    file = open(species_file,'r')
    lines = file.readlines()
    species_dict={}
    for line in lines:
        species_key = line.split(' ')[0]
        species_val = line.split(' ')[1:]
        species_val[-1] = species_val[-1].split("\n")[0]
        species_string = species_val[0]
        if len(species_val) > 1:
            for species in species_val[1:]:
                species_string = species_string + '-' + species
        species_dict[species_key] = species_string
    return(species_dict)


def make_descendants_dict(genes_file):
    file = open(genes_file,'r')
    lines = file.readlines()
    descendants_dict={}
    for line in lines:
        if len(line.split(' ')) > 2:
            gene = line.split(' ')[1]
            descendant = line.split(' ')[2]
            if '\n' in descendant:
                descendant = descendant.split('\n')[0]
            if '_' in descendant:
                descendant = descendant.split('_',1)[1]
            descendants_dict[gene] = descendant
    return(descendants_dict)


def replace_adjacencies(descendents_dict,adjacencies,geneorder_tab):
    gene1s=[]
    gene2s=[]
    for index, row in adjacencies.iterrows():
        if "|" in row["gene1"]:
            gene1_key=descendents_dict[row["gene1"]]
            if "|" in gene1_key:
                gene1_key = find_terminal_descendent(gene1_key,descendents_dict)
        elif "_" in row["gene1"]:
            #print(row["gene1"])
            gene1_key = row["gene1"].split('_',1)[1]
        new_gene1 = geneorder_tab[geneorder_tab["original_id"] == gene1_key]["gene"].iloc[0]
        gene1s.append(new_gene1)
        if "|" in row["gene2"]:
            gene2_key=descendents_dict[row["gene2"]]
            if "|" in gene2_key:
                gene2_key = find_terminal_descendent(gene2_key,descendents_dict)
        elif "_" in row["gene2"]:
            gene2_key = row["gene2"].split('_',1)[1]
        new_gene2 = geneorder_tab[geneorder_tab["original_id"] == gene2_key]["gene"].iloc[0]
        gene2s.append(new_gene2)
    adjacencies.drop("gene1", axis = 1, inplace = True)
    adjacencies.drop("gene2", axis = 1, inplace = True)
    adjacencies["gene1"] = gene1s
    adjacencies["gene2"] = gene2s
    #print(adjacencies)
    return(adjacencies)


def find_terminal_descendent(gene, descendants_dict):
    while '|' in gene:
        gene = descendants_dict[gene]
    return(gene)


def find_probable_ends(adjacencies,species_dict,focal_species):
    fronts = []
    backs =[]
    #for species in focal_species:
    for sp_num in list(set(adjacencies["species"])):
        if species_dict[str(sp_num)] in focal_species:
            subtable = adjacencies[adjacencies["species"] == int(sp_num)]
            front_end = find_front_end(subtable)
            back_end = find_back_end(subtable)
            for end in front_end:
                fronts.append(end)
            for end in back_end:
                backs.append(end)
    for i in fronts:
        if i in backs:
            backs.remove(i)
    for i in backs:
        if i in fronts:
            fronts.remove(i)
    fronts = list(set(fronts))
    backs = list(set(backs))
    return(fronts, backs)


def find_back_end(subtable):
    first_pair = [subtable["gene1"].iloc[0],subtable["gene2"].iloc[0]]
    second_pair = [subtable["gene1"].iloc[1],subtable["gene2"].iloc[1]]
    front_end = [x for x in first_pair if x not in second_pair]
    return(front_end)


def find_front_end(subtable):
    first_pair = [subtable["gene1"].iloc[-1],subtable["gene2"].iloc[-1]]
    second_pair = [subtable["gene1"].iloc[-2],subtable["gene2"].iloc[-2]]
    back_end = [x for x in first_pair if x not in second_pair]
    return(back_end)


def make_sorted_ends_list(fronts, backs):
    ends = deepcopy(fronts)
    for back in backs:
        ends.append(back)
    ends_dict ={}
    for i in list(set(ends)):
        ends_dict[i] = ends.count(i)
    ends = list(set(ends))
    ends = sorted(ends, key=lambda x : ends_dict[x],reverse=True)
    return(ends, ends_dict)


def parse_adjacencies_info(subtable):
    adjacencies_list = []
    adjacencies_objects = []
    for index, row in subtable.iterrows():
        entry = [row['gene1'], row['gene2']]
        adjacencies_list.append(entry)
        adjacencies_objects.append(row['gene1'])
        adjacencies_objects.append(row['gene2'])
    return(adjacencies_list, adjacencies_objects)


def find_starters_breaks(adjacencies_objects):
    adj_obj_counts = {}
    for obj in list(set(adjacencies_objects)):
        adj_obj_counts[obj] = adjacencies_objects.count(obj)
    starters = [obj for obj in list(adjacencies_objects) if adj_obj_counts[obj] == 1]
    breaks = [obj for obj in list(adjacencies_objects) if adj_obj_counts[obj] > 2]
    return(starters, breaks)


def merge_adjacencies(adjacencies_list):
    adj_len = len(adjacencies_list)
    changed_list = 0
    while changed_list != adj_len:
        adj_len = deepcopy(len(adjacencies_list))       
        for ref_adjacency in adjacencies_list:
            for test_adjacency in adjacencies_list[adjacencies_list.index(ref_adjacency)+1:]:
                if ref_adjacency[-1] == test_adjacency[0]:
                    for adj in test_adjacency[1:]:
                        ref_adjacency.append(adj)
                    adjacencies_list.remove(test_adjacency)
        for ref_adjacency in adjacencies_list:
            for test_adjacency in adjacencies_list[adjacencies_list.index(ref_adjacency)+1:]:
                if all((x in ref_adjacency) for x in test_adjacency):
                    adjacencies_list.remove(test_adjacency)
        changed_list = len(adjacencies_list)
    return(adjacencies_list)


def len_listobject(i):
  return len(i)


def identify_resolve_flips(adjacencies_list):
    front_end = [x[0] for x in adjacencies_list]
    back_end = [x[-1] for x in adjacencies_list]
    ends = deepcopy(front_end)
    for x in back_end:
        ends.append(x)
    doubles = [ gene for gene in ends if ends.count(gene) == 2]
    doubles = list(set(doubles))
    for double in doubles:
        pair = [adj for adj in adjacencies_list if double in adj]
        pair.sort(key=len_listobject)
        shorter = pair[0]
        adjacencies_list.remove(shorter)
        for adj in adjacencies_list:
            if adj[0]  == double or adj[-1] == double:
                if adj[0] == double:
                    for gene in shorter[1:]:
                        adj.insert(0, gene)
                    break
                elif adj[-1] == double:
                    shorter.reverse()
                    for gene in shorter[1:]:
                        adj.append(gene)
    return(adjacencies_list)
            

def merge_same_sides(adjacencies_list, fronts, backs):
    front_end = [x[0] for x in adjacencies_list]
    back_end = [x[-1] for x in adjacencies_list]
    for end in front_end:
        end_list = [x for x in adjacencies_list if x[0] == end]
        if len(end_list) > 1:
            for entry in end_list:
                adjacencies_list.remove(entry)
            end_list.sort(key=len_listobject)
            for entry in end_list:
                if back_end.count(entry[-1]) > 1:
                    end_list.remove(entry)
                    end_list.append(entry)
            for entry in end_list:
                if entry[-1] in backs:
                    end_list.remove(entry)
                    end_list.append(entry)
            focal = end_list[0]
            for string in end_list[1:]:
                for gene in string[1:]:
                    focal.append(gene)
            adjacencies_list.append(focal)
    for end in back_end:
        end_list = [x for x in adjacencies_list if x[-1] == end]
        if len(end_list) > 1:
            for entry in end_list:
                adjacencies_list.remove(entry)
            end_list.sort(key=len_listobject)
            for entry in end_list:
                if front_end.count(entry[0]) > 1:
                    end_list.remove(entry)
                    end_list.append(entry)
            for entry in end_list:
                if entry[0] in fronts:
                    end_list.remove(entry)
                    end_list.append(entry)
            focal = end_list[0]
            for string in end_list[1:]:
                string.reverse()
                for gene in string[1:]:
                    focal.insert(0,gene)
            adjacencies_list.append(focal)
    return(adjacencies_list)


def merge_opposite_sides(adjacencies_list):
    adj_len = len(adjacencies_list)
    changed_list = 0
    while changed_list != adj_len:
        adj_len = deepcopy(len(adjacencies_list))       
        for ref_adjacency in adjacencies_list:
            for test_adjacency in adjacencies_list[adjacencies_list.index(ref_adjacency)+1:]:
                if ref_adjacency[-1] == test_adjacency[0]:
                    for adj in test_adjacency[1:]:
                        ref_adjacency.append(adj)
                    adjacencies_list.remove(test_adjacency)
                if ref_adjacency[0] == test_adjacency[-1]:
                    test_adjacency.reverse()
                    for adj in test_adjacency[1:]:
                        ref_adjacency.insert(0,adj)
                    adjacencies_list.remove(test_adjacency)
        changed_list = len(adjacencies_list)
    return(adjacencies_list)


def find_front(adjacencies_list,fronts):
    possible_front = []
    for string in adjacencies_list:
        if string[0] in fronts or string[-1] in fronts:
            possible_front.append(string)
    for string in possible_front:
        if string[-1] in fronts:
            string.reverse()
    return(possible_front)


def find_back(adjacencies_list,backs):
    possible_back = []
    for string in adjacencies_list:
        if string[0] in backs or string[-1] in backs:
            possible_back.append(string)
    for string in possible_back:
        if string[0] in backs:
            string[0].reverse()
    return(possible_back)


def merge_flip_internal(adjacencies_list):
    adj_len = len(adjacencies_list)
    changed_list = 0
    while changed_list != adj_len:
        adj_len = deepcopy(len(adjacencies_list))       
        for ref_adjacency in adjacencies_list:
            for test_adjacency in adjacencies_list[adjacencies_list.index(ref_adjacency)+1:]:
                if (ref_adjacency[0] in test_adjacency or ref_adjacency[-1] in test_adjacency) or (test_adjacency[0] in ref_adjacency or test_adjacency[-1] in ref_adjacency):
                    if ref_adjacency[0] in test_adjacency or ref_adjacency[-1] in test_adjacency:
                        interior_match = test_adjacency
                        end_match = ref_adjacency
                    elif test_adjacency[0] in ref_adjacency or test_adjacency[-1] in ref_adjacency:
                        interior_match = ref_adjacency
                        end_match = test_adjacency
                    interior_index = deepcopy(adjacencies_list.index(interior_match))
                    adjacencies_list.remove(interior_match)
                    adjacencies_list.remove(end_match)
                    if end_match[0] in interior_match:
                        end_position = 'front'
                        interior_match_position = interior_match.index(end_match[0])
                        perserved_string = interior_match[0:interior_match_position]
                        flip_string = interior_match[interior_match_position:]
                        flip_string.reverse()
                        for gene in flip_string:
                            perserved_string.append(gene)
                        for gene in end_match[1:]:
                            perserved_string.append(gene)
                    elif end_match[-1] in interior_match:
                        end_position = 'back'
                        interior_match_position = interior_match.index(end_match[-1])
                        perserved_string = interior_match[interior_match_position+1:]
                        flip_string = interior_match[0:interior_match_position+1]
                        for gene in flip_string:
                            perserved_string.insert(0,gene)
                        end_match.reverse()
                        for gene in end_match[1:]:
                            perserved_string.insert(0,gene)
                    adjacencies_list.insert(interior_index,perserved_string)
        changed_list = len(adjacencies_list)
        return(adjacencies_list)
                    
                    
def filter_existing(adjacencies_list1):
    for merged_adj in adjacencies_list1:
        counter = 0
        for gene in merged_adj:
            for other_string in adjacencies_list1:
                if merged_adj != other_string and gene in other_string:
                    counter +=1
                    break
        if counter == len(merged_adj):
            adjacencies_list1.remove(merged_adj)
    return(adjacencies_list1)





def make_gene_strings(adjacencies_list,adjacencies,species_dict,focal_species):
    print("making gene string")
    fronts, backs = find_probable_ends(adjacencies,species_dict,focal_species)
    adjacencies_list1 = deepcopy(adjacencies_list)
    adjacencies_list1 = merge_adjacencies(adjacencies_list1)
    adjacencies_list1 = merge_opposite_sides(adjacencies_list1)
    adjacencies_list1 = identify_resolve_flips(adjacencies_list1)
    adjacencies_list1 = merge_same_sides(adjacencies_list1, fronts, backs)
    adjacencies_list1 = merge_opposite_sides(adjacencies_list1)
    adjacencies_list1 = filter_existing(adjacencies_list1)
    adjacencies_list1 = merge_flip_internal(adjacencies_list1)
    print(adjacencies_list1)
    if len(adjacencies_list1) != 1:
        front = find_front(adjacencies_list1,fronts)
        back = find_back(adjacencies_list1,backs)
        print(len(front),len(back))
        if len(front) == 1 and len(back) == 1 and front[0] != back[0]:
            final = front[0]
            for i in back[0]:
                final.append(i)
        elif len(front) == 1 and len(back) == 1 and front[0] == back[0]:
            final = front[0]
        elif len(front) == 1 and len(back) == 0 :
            final = front[0]
        elif len(front) == 0 and len(back) == 1 :
            final = back[0]
        else:
            final = []
    else:
        final = adjacencies_list1[0]
    print("hypothetically done making gene string")
    print(final)
    return(final)


def build_adjacencies_list(adjacencies,species_dict,focal_species, sp_num):
    subtable = adjacencies[adjacencies["species"] == sp_num]
    adjacencies_list, adjacencies_objects = parse_adjacencies_info(subtable)
    gene_strings = make_gene_strings(adjacencies_list,adjacencies,species_dict,focal_species)
    return(gene_strings)


def build_ancestors_table(geneorder_tab, adjacencies, species_dict, focal_species):
    final_table = pd.DataFrame(columns = ["molecules", "species", "gene", "start", "end", "strand", "orientation"])
    for sp_num in list(set(adjacencies["species"])):
        if species_dict[str(sp_num)] not in focal_species:
            print(sp_num)
            gene_strings = build_adjacencies_list(adjacencies,species_dict,focal_species, sp_num)
            #print(gene_strings)
            species_name = species_dict[str(sp_num)]
            species_list =  [species_name] * len(gene_strings)
            orthogroup = gene_strings
            strand = []
            for og in orthogroup:
                strand.append(list(geneorder_tab[geneorder_tab["gene"] == og]["strand"])[0])
            orientation = []
            for og in orthogroup:
                orientation.append(list(geneorder_tab[geneorder_tab["gene"] == og]["orientation"])[0])
            molecule = []
            for og in orthogroup:
                molecule.append(list(geneorder_tab[geneorder_tab["gene"] == og]["molecule"])[0])
            uniq_mol = list(set(molecule))
            for i in range(0, len(uniq_mol)):
                new_mol = species_name + "_mol" + str(i)
                molecule = [ new_mol if item == uniq_mol[i] else item for item in molecule]
            start = range(1,(len(gene_strings)*3),3)
            end = range(2,((len(gene_strings))*3),3)
            print(len(molecule), len(species_list), len(gene_strings), len(start), len(end), len(strand), len(orientation))
            final_subset = pd.DataFrame({"molecules":molecule,"species":species_list,"gene":gene_strings,"start":start,"end":end,"strand":strand,"orientation":orientation})
            final_table = pd.concat([final_table,final_subset])
    return(final_table)


def reformat_geneorder_tab(geneorder_tab):
    geneorder_tab = geneorder_tab[["molecule","species","gene","start","end","strand","orientation"]]
    newtable = pd.DataFrame(columns = ["molecules", "species", "gene", "start", "end", "strand", "orientation"])
    for i in list(set(geneorder_tab["species"])):
        subtable = deepcopy(geneorder_tab[geneorder_tab["species"] == i])
        subtable = subtable.sort_values(by=["molecule", "start"])
        start = list(range(1,(len(subtable)*3),3))
        end = list(range(2,((len(subtable))*3),3))
        subtable.drop("start", axis = 1, inplace = True)
        subtable.drop("end", axis = 1, inplace = True)
        final_subset = pd.DataFrame({"molecules":subtable["molecule"].to_list(),"species":subtable["species"].to_list(),"gene":subtable["gene"].to_list(),"start":start,"end":end,"strand":subtable["strand"].to_list(),"orientation":subtable["orientation"].to_list()})
        print(final_subset)
        newtable = pd.concat([newtable,final_subset])
    return(newtable)





             

                





########################################
################# CODE #################
########################################

species_file = decostar_dir + "/" + orthogroup + "/output/species.txt"
species_dict = read_decostar_species(species_file)

genes_file = decostar_dir + "/" + orthogroup + "/output/genes.txt"
genes_table = pd.read_csv(genes_file, ' ', names = ["species", "gene", "descendant1", "descendant2"])

geneorder_file = geneorder_dir + "/" + orthogroup + "/" + orthogroup + "_positions.csv"
geneorder_tab = pd.read_csv(geneorder_file,',')

adjacencies_file = decostar_dir + "/" + orthogroup + "/output/adjacencies.txt"
adjacencies = pd.read_csv(adjacencies_file, ' ', names = ["species", "gene1","gene2","number1","number2"])

focal_species = list(set(geneorder_tab["species"]))

descendants_dict = make_descendants_dict(genes_file)

adjacencies = replace_adjacencies(descendants_dict,adjacencies,geneorder_tab)

ancestor_table = build_ancestors_table(geneorder_tab, adjacencies, species_dict, focal_species)

existing_table = reformat_geneorder_tab(geneorder_tab)

final_table = pd.concat([ancestor_table,existing_table])


output_table_file = decostar_dir + "/" + orthogroup + "/output/ancestral_positions_table.csv"
final_table.to_csv(output_table_file, index=False)






