#!/usr/bin/env python

#import argparse
#import csv
from Bio.Seq import Seq
from Bio.Seq import IUPACData
import random
import numpy as np
from collections import Counter
import pandas as pd
from scipy.spatial import distance
import matplotlib.pyplot as plt
#from Bio import SeqIO


parser = argparse.ArgumentParser(description='Takes a fasta alignment and mapping file of individuals\' species assignments and will write a generax mapping file')
parser.add_argument("-f","--fasta",
					type=str,
					help="Fasta alignment")
parser.add_argument("-s","--starting_tree",
					type=str,
					default='',
					help="Starting tree")
parser.add_argument("-m","--model",
					type=str,
					help="Model for genetree estimation.")
parser.add_argument("-mp","--map",
					type=str,
					default='',
					help="Map file for genes and species.")
parser.add_argument("-o","--output_file",
					type=str,
					default='output',
					help="name of outputfile")
args=parser.parse_args()

fasta = args.fasta
starting_tree = args.starting_tree
model = args.model
map = args.map
output = args.output_file


# Write nexus file
outFile = open(output, "w")
outFile.write('[FAMILIES]\n')
outFile.write('-family_1\n')
if starting_tree != '':
	outFile.write('starting_gene_tree = ' + starting_tree + '\n')
outFile.write('alignment = ' + fasta + '\n')
if map != '':
	outFile.write('mapping = ' + map + '\n')
outFile.write('subst_model = ' + model + '\n')


outFile.close

#####################################################################
#########                    Functions                      #########
#####################################################################

## Iniital function to generate sequences
def generate_AA_seq(len):
    ## Define alphabet
	AA_aphabet = list(IUPACData.protein_letters)
    ## Select values from alphabet at random
	AAseq = ''.join(random.choice(AA_aphabet) for i in range(len-1))
    ## Make sure they start with M
	AAseq = ''.join(['M',AAseq])
    ## Return object
	AAseq = Seq(AAseq)
	return(AAseq)


## Differentiate sequences
## Input will be:
## 1) amino acid sequence
## 2) percent of differing sites
def diverge(AAseq, perc_div):
    
    ## Verify perc div is between 0 and 1
    if not (0 <= perc_div <= 1):
        raise ValueError("Divergence must be between 0 and 1.")
    
    ## Convert percent div into number of divergent sites
    num_differences = round(perc_div * len(AAseq))
    ## subtract 1 so first site can always be M
    num_differences = num_differences - 1
    ## Select sites to modify
    differing_sites = random.sample(range(1,len(AAseq)), num_differences)
    
    AAseq_div = list(AAseq)
    for site in differing_sites:
        ## Identify original residue to exclude
        original_residue = AAseq_div[site]
        ## Exclude original
        possible_mutations = [aa for aa in IUPACData.protein_letters if aa != original_residue]
        ## Replace with a random different AA
        AAseq_div[site] = random.choice(possible_mutations)
    
    ## Create new object
    AAseq_div = Seq(''.join(AAseq_div))
    return(AAseq_div)


## Function to generate gene family
## Takes starting sequence, generates N duplicates, 
## with mean variation MV drawn from a normal distribution
## with variance S
def create_gene_family(AAseq, N, MV, S):
    ## Make list of percent divergence values
    div = np.random.normal(MV, S, size=N)
    ## Create list of diverged sequences
    gene_family = [diverge(AAseq, div[i]) for i in range(N)]
    return(gene_family) 


## Function to differentiate gene family
## Can expand on this to allow variation by distribution
def diverge_gene_family(gene_family, M, S):
    ## Check that the gene family is >1 (otherwise use diverge seq)
    if not len(gene_family) > 1 :
        raise ValueError("gene family is less than or equal to 1. Try diverge instead")
    ## Make list of percent divergence values
    div = np.random.normal(M, S, size=len(gene_family))
    ## diverge each sequence by div amount
    diverged_seqs = [diverge(seq, div[index]) for index,seq in enumerate(gene_family)]
    return(diverged_seqs)


## Function to split an AA sequence into kmers of size K
def extract_kmers(AAseq, K):
    kmers = [AAseq[i:i+K] for i in range(0,len(AAseq)-K+1)]
    return(kmers)


## Function to make kmer counts dictionary given
##  a list of sequnces and kmer size K
def make_kmers_dict(AA_seq_list, K):
    kmers_lists = [extract_kmers(seq, K) for seq in AA_seq_list]
    ## Make list of all kmers across sequences
    All_kmer_list = [str(kmer) for kmers in kmers_lists for kmer in kmers]
    ## Now create the dictionary of counts and values
    kmer_counts = dict(Counter(All_kmer_list))
    return(kmer_counts)


## Function to make a kmers table of intersecting kmer counts
def intersect_kmer_dictionaries(kmer_list_1, kmer_list_2):
    kmers_table = pd.DataFrame.from_dict([kmer_list_1, kmer_list_2]).T  # Transpose to get keys as rows
    kmers_table.columns = ['Species1', 'Species2']  # Rename columns
    kmers_table.index.name = 'Key'  # Set index name
    kmers_table.fillna(0, inplace=True)
    return(kmers_table)


## Function to calculate dissimilarity
def kmer_bray_curtis(kmer_table):
    ## Convert kmer table to two vectors
    Species1 = kmer_table.iloc[:, 0]
    Species2 = kmer_table.iloc[:, 1]
    # Calculate Bray-Curtis distance 
    bray_curtis_dissimilarity = distance.braycurtis(Species1, Species2)
    # Convert to similarity (1 - distance)
    return(bray_curtis_dissimilarity)


## Function to add gene duplication
## Takes two starting gene families
## generates N duplications. Which gene family the  duplication
## occurs in and which gene are random. After the new sequences is 
## generated, it is diverged from its original sequence by 
## a divergence percentage value drawn from a normal 
## distribution with mean (MD) and variance S
def gene_duplications(gene_family_1, gene_family_2, N, MD, S):
    new_family_1 = gene_family_1.copy()
    new_family_2 = gene_family_2.copy()
    counter=1
    while counter <= N:
        counter += 1
        family = random.choice([new_family_1, new_family_2])   
        gene = random.choice(family)
        if MD == 0:
            new_seq = gene  
        else:
            div = np.random.normal(MD, S)
            new_seq = diverge(gene, div)         
        family.append(new_seq)
    return(new_family_1, new_family_2)



## Function to add gene loss
## Takes two starting gene families
## generates N gene losses. Which gene family the loss
## occurs in and which gene are random.
def gene_loss(gene_family_1, gene_family_2, N):
    new_family_1 = gene_family_1.copy()
    new_family_2 = gene_family_2.copy()
    counter=1
    while counter <= N:
        counter += 1
        family = random.choice([new_family_1, new_family_2])   
        gene = random.choice(family)
        family.remove(gene)         
    return(new_family_1, new_family_2)


## Function to add gene loss
## Takes two starting gene families
## generates N gene losses. Which gene family the loss
## occurs in and which gene are random.
def gene_expression(gene_family_1, 
                    gene_family_2, 
                    expression_mean,
                    expression_var, 
                    difference_mean,
                    difference_var):
    new_family_1 = gene_family_1.copy()
    new_family_2 = gene_family_2.copy()
    for gene in list(set(new_family_1)):
        N_expression_gf1 = np.random.normal(expression_mean, 
                                        expression_var)
        for _ in range(int(N_expression_gf1-1)):
            new_family_1.append(gene)
    for gene in list(set(new_family_2)):
        N_expression_gf1 = np.random.normal(expression_mean, 
                                        expression_var)
        diff = np.random.normal(difference_mean,
                                difference_var)
        ## We'll just add so that we avoid negative numbers
        ## gene family 2 will always be more highly expressed
        ## but that should not matter since we aren't interested
        ## in directionality here.
        N_expression_gf2 = N_expression_gf1 + diff
        for _ in range(int(N_expression_gf2-1)):
            new_family_2.append(gene)        
    return(new_family_1, new_family_2)



#####################################################################
#########                       Code                       #########
#####################################################################

AA_seq_SVMPlen = generate_AA_seq(600)
AA_gene_family1 = create_gene_family(AA_seq_SVMPlen, 5, 0.10, 0.01)
AA_gene_family2 = diverge_gene_family(AA_gene_family1, 0.10, 0)


gene_family1_kmers = make_kmers_dict(AA_gene_family1,20)
gene_family2_kmers = make_kmers_dict(AA_gene_family2,20)

## 1 duplication, no divergence
gene_family1_mod, gene_family2_mod = gene_duplications(AA_gene_family1, AA_gene_family2, N=1, MD=0, S=0)

SVMPlen_kmers_table = intersect_kmer_dictionaries(gene_family1_kmers, gene_family2_kmers)

kmer_bray_curtis(SVMPlen_kmers_table)


###################################################################
## Tests 
## 1) Functional div across random AA variation
perc_div = []
functional_divergences =[]
## set range of divergence levels
for i in np.arange(0.01, 1.00, 0.01):
    counter = 0
    ## for each value of diverence do 10 replicates
    while counter < 10:
        counter += 1
        ## Set starting seq
        seq1 = AA_seq_SVMPlen
        ## Create second sequence diverged by i%
        seq2 = diverge(seq1, i)
        seq1 = [seq1]
        seq2 = [seq2]
        ## Pull kmers
        gene_family1_kmers = make_kmers_dict(seq1,20)
        gene_family2_kmers = make_kmers_dict(seq2,20)
        kmer_table = intersect_kmer_dictionaries(gene_family1_kmers, gene_family2_kmers)
        ## Calc functional divergence
        funct_div = kmer_bray_curtis(kmer_table)
        perc_div.append(i)
        functional_divergences.append(funct_div)

## Plot
plt.scatter(perc_div, functional_divergences, s=30)
plt.xlabel('Percent amino acid divergence')
plt.ylabel('Single sequence functional divergence')
plt.title('Functional divergence versus amino acid divergence')
plt.savefig('seq_aa_div.png')
plt.savefig('seq_aa_div.svg')
plt.clf()
#################################################################



## 2) Functional div across random AA variation with sequences of different lengths
perc_div = []
functional_divergences =[]
seq_size=[]
## set range of divergence levels
for size in [139, 259, 600]:
    for i in np.arange(0.01, 1.00, 0.01):
        counter = 0
        ## for each value of diverence do 10 replicates
        while counter < 10:
            counter += 1
            ## Set starting seq
            seq1 = generate_AA_seq(size)
            ## Create second sequence diverged by i%
            seq2 = diverge(seq1, i)
            seq1 = [seq1]
            seq2 = [seq2]
            ## Pull kmers
            gene_family1_kmers = make_kmers_dict(seq1,20)
            gene_family2_kmers = make_kmers_dict(seq2,20)
            kmer_table = intersect_kmer_dictionaries(gene_family1_kmers, gene_family2_kmers)
            ## Calc functional divergence
            funct_div = kmer_bray_curtis(kmer_table)
            perc_div.append(i)
            functional_divergences.append(funct_div)
            seq_size.append(size)


colors = {139: '#461763', 259: '#26ad81', 600: '#e9e223'}
## Plot
for length in [139, 259, 600]:
    idx = [i for i, g in enumerate(seq_size) if g == length]
    plt.scatter([perc_div[i] for i in idx], [functional_divergences[i] for i in idx],
                c=colors[length], label=length, s=30)

plt.xlim(0, 0.6)
plt.xlabel('Percent amino acid divergence')
plt.ylabel('Single sequence functional divergence')
plt.title('Functional divergence versus amino acid divergence')
plt.legend(title="AA seq size")
plt.savefig('seq_aa_len_div.png')
plt.savefig('seq_aa_len_div.svg')
plt.clf()



#################################################################


#################################################################
## 3) Functional div across gene duplication variation
perc_div = []
functional_divergences =[]
n_gene_duplicates = []
## set range of gene duplications
for ngene in range(0,21):
    ## set range of divergence levels
    ## here i will track divergence
    for i in np.arange(0.01, 0.40, 0.01):
        counter = 0
        while counter < 10:
            counter += 1
            seq1 = AA_seq_SVMPlen
            ## We will assume a flat 10% divergent sites.
            seq2 = diverge(seq1, 0.1)
            seq1 = [seq1]
            seq2 = [seq2]
            ## Here we will do ngene (N) duplications
            ## We will assume a constant rate of divergence among
            ## all genes, thus duplicate genes will have the same
            ## percent of divergent sites as the starting gene families (MD=i)
            gene_family1, gene_family2 = gene_duplications(seq1, seq2, N=ngene, MD=i, S=0)
            gene_family1_kmers = make_kmers_dict(gene_family1,20)
            gene_family2_kmers = make_kmers_dict(gene_family2,20)
            kmer_table = intersect_kmer_dictionaries(gene_family1_kmers, gene_family2_kmers)
            funct_div = kmer_bray_curtis(kmer_table)
            perc_div.append(i)
            functional_divergences.append(funct_div)
            n_gene_duplicates.append(ngene)


## Plot

#plt.scatter(perc_div, functional_divergences, s=30,  c=n_gene_duplicates, cmap='viridis')
sc = plt.scatter(n_gene_duplicates, functional_divergences, s=30,  c=perc_div, cmap='viridis')
cbar = plt.colorbar(sc, label='Color scale')  # Create colorbar
cbar.set_label('Proportion of differing sites') # Set color scale label
plt.xlim(-1, 21)
plt.xlabel('Number of gene duplicates')
plt.ylabel('Functional divergence')
plt.title('Functional divergence versus amino acid divergence with gene duplication',
          fontsize=10)
plt.savefig('seq_gene_aa_div.png')
plt.savefig('seq_gene_aa_div.svg')
plt.clf()


#################################################################


## 4) Test the effect of expression divergence
#perc_div = []
functional_divergences =[]
diff_in_expression = []
## set range of divergence levels
## here i will track divergence
for diff in [6,12,25,50,100,200,400,800,1600,3200,6400,
            12800,25600,51200,102400]:
    counter = 0
    while counter < 10:
        counter += 1
        seq1 = AA_seq_SVMPlen
        AA_seq_SVMPlen = generate_AA_seq(600)
        AA_gene_family1 = create_gene_family(AA_seq_SVMPlen, 5, 0.10, 0.01)
        AA_gene_family2 = diverge_gene_family(AA_gene_family1, 0.10, 0)
        ## We will use our gene families 1 and 2
        ## which are SVMP len sequences
        ## 5 genes in each family, each with 10%
        ## site variation among paralogs
        ## and then 10% site variation between the two families
        gene_family_1_expr, gene_family_2_expr = gene_expression(AA_gene_family1, 
                                                                AA_gene_family2,
                                                                100, 0, diff, 0)
        gene_family1_expr_kmers = make_kmers_dict(gene_family_1_expr,20)
        gene_family2_expr_kmers = make_kmers_dict(gene_family_2_expr,20)
        kmer_table = intersect_kmer_dictionaries(gene_family1_expr_kmers, gene_family2_expr_kmers)
        funct_div = kmer_bray_curtis(kmer_table)
        #perc_div.append(i)
        functional_divergences.append(funct_div)
        diff_in_expression.append(diff)


log_diff_in_expression = np.log10(diff_in_expression)
## Plot
#plt.scatter(perc_div, functional_divergences, s=30,  c=n_gene_duplicates, cmap='viridis')
plt.scatter(log_diff_in_expression, functional_divergences, s=30)
#cbar = plt.colorbar(sc, label='Color scale')  # Create colorbar
#cbar.set_label('Proportion of differing sites') # Set color scale label
#plt.xlim(3, 12)
plt.xlabel('Difference Between Gene Family Expression (log10(diff))')
plt.ylabel('Functional Divergence')
plt.title('Functional divergence versus difference in gene family expression',
          fontsize=10)
plt.savefig('expr_aa_div.png')
plt.savefig('expr_aa_div.svg')
plt.clf()



