#!/usr/bin/env python

# Additional software necessary to run this:

# (1) biopython
# (2) pandas
# (3) dfply
# (4) mafft
# (5) iqtree
# (6) generax v.XXX

import argparse
import copy
from re import sub
import threading
from ete3 import Tree
from ete3 import PhyloTree
from itertools import combinations
import subprocess as sp
import pandas as pd
import csv
#from ete3 import Phyloxml, phyloxml

from ctypes import alignment
import sys, os, shutil
import datetime as dt
#import numpy as np
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
parser.add_argument("-st","--species_tree",
					type=str,
                    default="0_Inputs/species_tree.nwk",
					help="Species tree file for orthofinder. Newick format")
parser.add_argument("-rd","--reconciled_dir",
					type=str,
                    default="3.5_OG_reconciliations",
					help="Directory containing generax input and output files")
parser.add_argument("-md","--map_dir",
					type=str,
                    default="3.1_Species_maps",
					help="Directory containing sequence")
parser.add_argument("-gd","--genomes_dir",
					type=str,
                    default="0_Inputs/genomes",
					help="Directory containing the genomes for focal taxa. Strictly used to identify focal taxa")
parser.add_argument("-od","--output_dir",
					type=str,
                    default="3.6_OG_classifications",
					help="Directory of where to write a folder for this orthogroup and subsequent output files")
args = parser.parse_args()

########################################
################# SETUP ################
########################################

orthogroup = args.orthogroup
reconciled_dir = args.reconciled_dir
map_dir = args.map_dir
species_tree = args.species_tree
genomes_dir = args.genomes_dir
output_dir = args.output_dir
#out = args.output_directory


########################################
############## FUNCTIONS ###############
########################################

def find_species(genomes_dir) :
	species_list = list(os.listdir(genomes_dir))
	species_list = [x.split('.')[0] for x in species_list]
	return(species_list)


def parse_sp_name(node_name):
    return node_name.split("_")[0]


def read_and_reconcile_tree(orthogroup, reconciled_dir, map):
    tree_string = reconciled_dir + "/" + orthogroup + "/" + orthogroup + "_output/results/family_1/geneTree.newick"
    tree = PhyloTree(tree_string)
    ## Renames tree leaves to have Species
    for leaf in tree.get_tree_root():
        new_name = map[map["Sequence"] == leaf.name]["New_name"].iloc[0]
        leaf.name = new_name
    tree.set_species_naming_function(parse_sp_name)
    #genetree = PhyloTree(tree)
    #st = PhyloTree(st)
    recon_tree, events = tree.reconcile(st)
    return(tree)


def modify_classifications(tree, st, focal_species):
    tree = modify_in_paralog(tree)
    identify_and_resolve_topology_problems(tree,st)
    identify_putative_false_orthology(tree)
    final_paralogs_list,nhx_tree = fix_in_paralogs_define_groups(tree, focal_species)
    final_paralogs_list = filter_nonfocal_species(final_paralogs_list, focal_species)
    return(final_paralogs_list, nhx_tree)


def identify_putative_false_orthology(tree):
    for node in tree.traverse():
        if (hasattr(node,'evoltype') and node.evoltype == 'D'):
            for child in node.get_children():
                if hasattr(child,'evoltype') and child.evoltype == 'S' and len(list(child.get_species())) != 1:
                    ref = child.dist
                    print(ref)
                    ## get_children means we are only looking at two nodes. Will also speed things up.
                    for descendent in child.get_children():
                        if long_br(descendent,ref) and descendent.up.evoltype != 'SD':
                            print('')
                            print("ref ", ref)
                            print("made a change based on ")
                            print(descendent)
                            descendent.up.evoltype = 'D'
                            print(descendent.up)
                            #print('')


def long_br(descendent,ref):
    if ref == 0.0:
        ref = 0.000001
    test = descendent.dist
    if test == 0.0:
        test = 0.000001
    print("test", test)
    val = test /ref
    print("val", val)
    if val > 5 and ref != 0.000001:
        return True
    else:
        return False


def identify_and_resolve_topology_problems(tree,st):
    for node in tree.traverse():
        if hasattr(node, 'evoltype') and node.evoltype == 'D' :
            for child in node.get_children():
                if hasattr(child,'evoltype') and child.evoltype == 'S' :
                    #print("running test")
                    if short_br(node,child) and species_rep_test(node):
                        #node_to_fix = node.detach()
                        #print(node)
                        fix_node(tree,node,st)


def short_br(node,child):
    ref = node.dist
    if ref == 0.0:
        ref = 0.00000000001
    test = node.get_distance(child)
    if test == 0.0:
        test = test = 0.00000000001
    val = test /ref
    #print("parent is ", str(ref))
    #print("child is ",str(test))
    #print("test ratio ", str(val))
    if val < 0.1:
        return True
    else:
        return False


## Tests number of sepcies representatives below a given node accounting for in paralogs.
## Passes (True) if only 1 per species. Fails (False) if more than 1 after account for in-paralogs.
def species_rep_test(node):
    species = [x.species for x in node.get_leaves()]
    sd_species = [x.species for x in node.traverse() if hasattr(node,"evoltype") and node.evoltype == "SD"]
    if len(sd_species) != 0:
        for i in sd_species:
            species.remove(i)
    if any([species.count(x) > 1 for x in species]):
        return False
    else:
        return True


def fix_node(tree,bad_node,st):
    node = copy.deepcopy(bad_node)
    spare_st = copy.deepcopy(st)
    represented_species = [x for x in spare_st.get_leaf_names() if x in node.get_species()]
    spare_st.prune(represented_species)
    outgroup_species = find_reference_outgroup(spare_st)
    species_instances = [x.species for x in node.get_leaves()]
    if len(outgroup_species) > 1 or species_instances.count(outgroup_species[0]) > 1:
        species_leaves = search_for_species(node,outgroup_species)
        out_anc = node.get_common_ancestor(species_leaves)
        node.set_outgroup(out_anc)
    else:
        out_leaf = search_for_species(node,outgroup_species)[0]
        node.set_outgroup(out_leaf)
    rec_node, node_events = node.reconcile(spare_st)
    #print(node)
    #for i in node.traverse():
        #if hasattr(i,"evoltype"):
            #print(i.evoltype)
    bad_node.up.add_child(node)
    bad_node.detach()
    

def search_for_species(tree, species):
    node_list = []
    for entry in species:
        for node in tree.traverse():
            if hasattr(node,"species") and node.species == entry:
                node_list.append(node)
    return(node_list)


def find_reference_outgroup(st):
    outgroup = list(st.get_leaf_names())
    for child in st.get_children():
        if child.is_leaf() == True:
            outgroup = [child.name]
            break
        elif len(list(child.get_leaf_names())) < len(outgroup):
            outgroup = list(child.get_leaf_names())
    return(outgroup)


def modify_in_paralog(tree):
    for node in tree.traverse():
        if hasattr(node, "evoltype") and node.evoltype == "D":
            if len(list(node.get_species())) == 1:
                node.evoltype = "SD"
                sd_species = list(node.get_species())[0]
                node.species = sd_species
    return(tree)



def fix_in_paralogs_define_groups(tree,focal_species):
    #print(tree)
    nhx_tree = copy.deepcopy(tree)
    for node in nhx_tree.traverse():
        if hasattr(node, "evoltype") and node.evoltype == "SD":
            node.evoltype = "D"
    final_paralog_list = []
    for subtree in tree.split_by_dups(autodetect_duplications=False):
        SD_events = count_SD_events(subtree)
        print(SD_events)
        print(subtree)
        if SD_events == 0:
            add_evolevent_feature(subtree,focal_species)
            final_paralog_list.append(subtree)
        else:
            in_paralogs = SD_graft(subtree)
            for paralog in in_paralogs:
                final_paralog_list.append(paralog)
    return(final_paralog_list , nhx_tree)


def count_SD_events(subtree):
    counter = 0
    for node in subtree.traverse():
        if hasattr(node, "evoltype") and node.evoltype == "SD":
            counter += 1
    return(counter)
        

def SD_graft(subtree):
    node_list = []
    if len(subtree.get_species()) == 1:
        for node in subtree.traverse():
            if hasattr(node, "evoltype") and node.evoltype == "SD":
                node.evoltype = "D"
        ##print("test")
        paralogs = subtree.detach()
        ##print(paralogs)
        for paralog in paralogs.split_by_dups(autodetect_duplications=False):
            ##print(9)
            ##print(paralog) 
            species = [x.species for x in paralog.get_leaves()]
            paralog.add_feature("evolevent","duplication")
            paralog.add_feature("duplicated_taxon",species)
            node_list.append(paralog)
    else:
        nnodespecies = len(subtree.get_species())
        while len([x.species for x in subtree.get_leaves()]) != nnodespecies:
            for node in subtree.traverse():
                ##print(1)
                #while hasattr(node, "evoltype") and node.evoltype == "SD" :
                if hasattr(node, "evoltype") and node.evoltype == "SD" :
                    nnodespecies = len(list(set([x.species for x in node.get_leaves()])))
                    farthest, test = node.get_farthest_leaf()
                    for leaf in node.get_leaves():
                        ##print(2)
                        dist = node.get_distance(leaf)
                        ##print(dist)
                        ##print(test)
                        if dist <= test:
                            dist = test
                            chosen_leaf = leaf
                            chosen_leaf_copy = copy.deepcopy(chosen_leaf)
                            chosen_leaf_copy.dist = dist
                        ## now we have the leaf to graft to node
                    print(3)
                    node.up.add_child(chosen_leaf_copy)
                    ### Here is our break.
                    paralogs = node.detach()
                    nspecies = [x.species for x in subtree.get_leaves()]
                    print(4)
                    if all([nspecies.count(x) == 1 for x in focal_species]):
                        print(5)
                        subtree.add_feature("evolevent","conserved")
                    elif any([nspecies.count(x) < 1 for x in focal_species]):
                        print(6)
                        lost_species_list = [x for x in focal_species if nspecies.count(x) < 1]
                        subtree.add_feature("evolevent","loss")
                        subtree.add_feature("lost_taxa",lost_species_list)
                    print(7)
                    break
            for node in paralogs.traverse():
                print(8)
                if hasattr(node, "evoltype") :
                    node.evoltype ="D"
            print(paralogs.evoltype)
            for paralog in paralogs.split_by_dups(autodetect_duplications=False):
                print(9)
                print(paralog.get_leaf_names())
                print(list(chosen_leaf_copy.get_leaf_names())[0])
                if list(paralog.get_leaf_names())[0] != list(chosen_leaf_copy.get_leaf_names())[0]:
                    species = [x.species for x in paralog.get_leaves()]
                    paralog.add_feature("evolevent","duplication")
                    paralog.add_feature("duplicated_taxon",species)
                    node_list.append(paralog)
                print(10)
            nnodespecies = len(subtree.get_species())
        node_list.append(subtree)
        print(subtree)
    return(node_list)


def add_evolevent_feature(node,focal_species):
    #allspecies = list(tree.get_species())
    nspecies = [x.species for x in node.get_leaves()]
    if all([nspecies.count(x) == 1 for x in focal_species]):
        node.add_feature("evolevent","conserved")
    elif any([nspecies.count(x) < 1 for x in focal_species]):
        lost_species_list = [x for x in focal_species if nspecies.count(x) < 1]
        node.add_feature("evolevent","loss")
        node.add_feature("lost_taxa",lost_species_list)
    elif any([nspecies.count(x) > 1 for x in focal_species]):
        for child in node.traverse():
            if hasattr(child,"evolevent") and child.evolevent == "SD":
                child.add_feature("evolevent","duplication")
    else:
        print("uhhh")




def filter_nonfocal_species(subtree_list, focal_species):
    for subtree in subtree_list:
        species = list(subtree.get_species())
        if not any(x in species for x in focal_species):
            subtree_list.remove(subtree)
        elif not all(x in species for x in focal_species):
            keep_leaves = []
            for node in subtree:
                if hasattr(node, "species") and node.species in focal_species:
                    keep_leaves.append(node.name)
            subtree.prune(keep_leaves)
    return(subtree_list)
                    

def build_orthologs_output(paralogs_list):
    orthologs_out = []
    counter = 0
    for paralog in paralogs_list:
        counter += 1
        gene_name = orthogroup + "-Gene-" + str(counter)
        string = ''
        for leaf_name in paralog.get_leaf_names():
            string = string + leaf_name
            if leaf_name != paralog.get_leaf_names()[-1]:
                string = string + ' '
        orthologs_out.append([gene_name,string])
    return(orthologs_out)


def build_ortholog_keys(paralogs_list):
    ortholog_keys = []
    counter = 0
    for paralog in paralogs_list:
        counter += 1
        gene_name = orthogroup + "-Gene-" + str(counter)
        for leaf_name in paralog.get_leaf_names():
            entry = [leaf_name,gene_name]
            ortholog_keys.append(entry)
    return(ortholog_keys)


def build_classes_table(paralogs_list, focal_species):
    classes = []
    counter = 0
    species_comparisons = [comb for comb in combinations(focal_species, 2)]
    header = ['paralog']
    for comp in species_comparisons:
        comp_string = comp[0] + '-' + comp[1]
        header.append(comp_string)
    classes.append(header)    
    for paralog in paralogs_list:
        counter += 1
        gene_name = orthogroup + "-Gene-" + str(counter)
        entry = [gene_name]
        if paralog.evolevent == "conserved":
            n_conserved = ["conserved"] * len(species_comparisons)
            for i in n_conserved:
                entry.append(i)
        elif paralog.evolevent == "loss":
            for comp in species_comparisons:
                if comp[0] not in paralog.lost_taxa and comp[1] not in paralog.lost_taxa:
                    entry.append("conserved")
                elif comp[0] in paralog.lost_taxa and comp[1] in paralog.lost_taxa:
                    entry.append("NA")
                else:
                    entry.append("loss")
        elif paralog.evolevent == "duplication":
            for comp in species_comparisons:
                if comp[0] in paralog.duplicated_taxon or comp[1] in paralog.duplicated_taxon:
                    entry.append("duplication")
                else:
                    entry.append("NA")
        classes.append(entry)
    return(classes)


def fix_nhx_names(nhx_tree,map):
    for leaf in nhx_tree.get_leaves():
        new_name = map[map["New_name"] == leaf.name]["Sequence"].iloc[0]
        leaf.name = new_name
    return(nhx_tree)


def fix_paralog_names(paralog_list, map):
    for paralog in paralog_list:
        for leaf in paralog.get_leaves():
            new_name = map[map["New_name"] == leaf.name]["Sequence"].iloc[0]
            leaf.name = new_name
    return(paralog_list)


def modify_nhx(nhx,focal_species):
    new_tree = copy.deepcopy(nhx)
    keep_leaves = [leaf.name for leaf in new_tree.get_leaves() if leaf.species in focal_species]
    new_tree.prune(keep_leaves)
    return(new_tree)


def decostar_sp_tree(species_tree,focal_species):
    keep_leaves = [leaf.name for leaf in species_tree.get_leaves() if leaf.name in focal_species]
    edited_species_tree = copy.deepcopy(species_tree)
    edited_species_tree.prune(keep_leaves)
    return(edited_species_tree)


def add_anc_names(st):
    possible_species = st.get_leaf_names()
    counter = 0
    for node in st.traverse():
        if node.species == '':
            sp_name = "anc" + str(counter)
            node.add_feature("S",sp_name)
            counter += 1
        else:
            putative_species = [x for x in possible_species if node.species in x]
            if len(putative_species) == 1:
                node.add_feature("S", putative_species[0])
    return(st)



#def make_anc_dict(focal_st):
#    anc_dict = {}
#    for node in focal_st.traverse():
#        species = tuple(node.get_leaf_names())
#        anc_dict[species] = node.S
#    return(anc_dict)



def make_decostar_nhx(nwk_tree, focal_st):
    decostar_nhx = copy.deepcopy(nwk_tree)
    for node in decostar_nhx.traverse():
        if hasattr(node, 'evoltype') and node.evoltype == "D":
            node.add_feature("D",'Y')
        else:
            node.add_feature("D","N")
    for node in decostar_nhx.traverse():
        all_observed_species = [x.split('_')[0] for x in node.get_leaf_names()]
        uniq_observed_species = list(set(all_observed_species))
        if node.is_leaf():
            anc = node.get_leaf_names()[0].split('_')[0]
            node.add_feature("S",anc)
            node.add_feature("geneName",node.get_leaf_names()[0])
        elif node.evoltype == "S":
            anc = focal_st.get_common_ancestor(uniq_observed_species).S
            node.add_feature("S",anc)
        elif node.evoltype == "D" and len(uniq_observed_species) == 1:
            anc = uniq_observed_species[0]
            node.add_feature("S",anc)
        elif node.evoltype == "D" :
            anc = focal_st.get_common_ancestor(uniq_observed_species).S
            node.add_feature("S",anc)
    return(decostar_nhx)
    


########################################
################# CODE #################
########################################

st = PhyloTree(species_tree)

focal_species = find_species(genomes_dir)
if '' in focal_species:
    focal_species.remove('')
if '.dir' in focal_species:
    focal_species.remove('.dir')


focal_st = copy.deepcopy(st)
focal_st.prune(focal_species)
focal_st = add_anc_names(focal_st)

map_file = map_dir + "/" + orthogroup + '.map'
map = pd.read_csv(map_file,",")

tree = read_and_reconcile_tree(orthogroup, reconciled_dir, map)

paralogs_list,nhx_tree = modify_classifications(tree, st, focal_species)

nwk_tree = modify_nhx(nhx_tree,focal_species)

decostar_nhx = make_decostar_nhx(nwk_tree, focal_st)
print(decostar_nhx.is_root())


nhx_tree = fix_nhx_names(nhx_tree, map)

paralogs_list = fix_paralog_names(paralogs_list, map)

orthologs_table = build_orthologs_output(paralogs_list)

ortholog_keys = build_ortholog_keys(paralogs_list)

classes_table = build_classes_table(paralogs_list, focal_species)

edited_species_tree = decostar_sp_tree(st,focal_species)

########################################
############### Outputs ################
########################################

mk_output = "mkdir " + output_dir + "/" + orthogroup
sp.call(mk_output, shell=True)


out_tree = output_dir + "/" + orthogroup + "/" + orthogroup + ".nhx"
nhx_tree.write(features = ["species","evoltype"], outfile=out_tree, format_root_node=True)

out_nwk = output_dir + "/" + orthogroup + "/" + orthogroup + "_decostar.nwk"
nwk_tree.write(format=9, outfile=out_nwk, format_root_node=True)
out_codeml = output_dir + "/" + orthogroup + "/" + orthogroup + "_codeml.nwk"
nwk_tree.write(outfile=out_codeml, format_root_node=True)


out_deco = output_dir + "/" + orthogroup + "/" + orthogroup + "_decostar.nhx"
decostar_nhx.write(format=9,features = ["S","D"], outfile=out_deco, format_root_node=True)

out_expr = output_dir + "/" + orthogroup + "/" + orthogroup + "_expr.nhx"
decostar_nhx.write(format=0,features = ["S","D"], outfile=out_expr, format_root_node=True)


out_species = output_dir + "/" + orthogroup + "/" + orthogroup + "_pruned_species.nwk"
edited_species_tree.write(format=9,outfile=out_species, format_root_node=True)


orthologs_outfile = output_dir + '/' + orthogroup + "/" + orthogroup + '_orthogroup.csv'
with open(orthologs_outfile, 'w') as csv_file:
    csv_writer = csv.writer(csv_file, delimiter = ',')
    for row in orthologs_table:
        csv_writer.writerow(row)		
csv_file.close()


ortho_keys_outfile = output_dir + '/' + orthogroup + "/" + orthogroup + '_genes.csv'
with open(ortho_keys_outfile, 'w') as csv_file:
    csv_writer = csv.writer(csv_file, delimiter = ',')
    csv_writer.writerow(['gene','orthogroup'])
    for row in ortholog_keys:
        csv_writer.writerow(row)		
csv_file.close()


classes_table_outfile = output_dir + '/' + orthogroup + "/" + orthogroup + '_classes.csv'
with open(classes_table_outfile, 'w') as csv_file:
    csv_writer = csv.writer(csv_file, delimiter = ',')
    for row in classes_table:
        csv_writer.writerow(row)		
csv_file.close()





#tree_string = reconciled_dir + "/" + orthogroup + "/" + orthogroup + "_output/results/family_1/geneTree.newick"
#    tree = PhyloTree(tree_string)