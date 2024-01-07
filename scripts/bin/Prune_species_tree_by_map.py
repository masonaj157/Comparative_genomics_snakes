#!/usr/bin/env python

import argparse
import csv
from ete3 import Tree

parser = argparse.ArgumentParser(description='Takes a fasta alignment and mapping file of individuals\' species assignments and will write a generax mapping file')
parser.add_argument("-st","--species_tree",
					type=str,
					help="Species tree to translate. Script assumes your tree has leaf names, branch lengths, and it supports branch supports though they are not required")
parser.add_argument("-m","--map",
					type=str,
					help="Map file formatted for generax (Species:Indiv;Indiv). Assumes that species not represented in the gene tree are not included in the map")
parser.add_argument("-o","--output",
					type=str,
					default='output',
					help="name of output tree file written in newick format")
args=parser.parse_args()

tree = args.species_tree
map = args.map
output = args.output

#Parse tree

tree = Tree(tree)

#Parse map
#As noted above, we are assuming species not represented in the gene tree are not represented in the map
species = []
file = open(map, "r")
file = file.readlines()
for line in file:
	line = line.split(':')[0]
	species.append(line)


tree.prune(species, preserve_branch_length=True)


tree.write(format=0, outfile=output)

		