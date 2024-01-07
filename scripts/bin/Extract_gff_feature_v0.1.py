#!/usr/bin/env python


import argparse
import csv
from Bio import SeqIO
from Bio import Seq
from Bio.SeqRecord import SeqRecord

parser = argparse.ArgumentParser(description='A script to extract gene or CDS sequences from a genome given a genome fasta and gff. Sadly quite slow, but I\'m no computational biologist. ')
parser.add_argument("-s","--seqfile",
					type=str,
					default='',
					help="Sequence file in fasta format")
parser.add_argument("-g","--gff",
					type=str,
					default='',
					help="gff3 file containing genes as a type")
parser.add_argument("-f","--feature",
					type=str,
					help="feature you would like to collect in a fasta file. Must be a string matched in your gff")
parser.add_argument("-o","--output",
					type=str,
					default='',
					help="name of the output CDS file")
parser.add_argument("-r","--remove_cds_label",
					type=str,
					default=True,
					help="Should the \':cds\' string be removed? True or False. Default is true to correspond to featureCounts output.")
args=parser.parse_args()

seqfile = args.seqfile
gff = args.gff
feature = args.feature
output = args.output
remove_cds_label = args.remove_cds_label

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
        self.attributes = parse_GFF_attributes(gff_entry[8])
        if self.attributes.ID == []:
            self.attributes.ID = ['unnamed']


class gff_attributes:
    def __init__(self, gff_attributes):
        self.ID = gff_attributes[0]
        self.Name = gff_attributes[1]
        self.Alias = gff_attributes[2]
        self.Note = gff_attributes[3]



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
    attributes_list = gff_attributes([ID, Name, Alias, Note])
    return(attributes_list)



def pull_gene(sequence_list, gff_list, feature) :
    gff_list = [ entry for entry in gff_list if entry.type == 'gene']
    feature_gffs = [entry for entry in gff_list if entry.type == feature]
    feature_seqs = []
    for entry in feature_gffs:
        new_seq = find_gene_seq(sequence_list,entry)
        if not isinstance(new_seq, str) :
            feature_seqs.append(new_seq)
    return(feature_seqs)


def test_for_multiple_gff_matches(seq, entry):
    if len(seq) > 1:
        print('so, somehow there is a gff that matches more than one sequence')
        print(entry.attributes.ID)
        return True
    else:
        return False
        

def find_gene_seq(sequence_list,entry):
    if entry.attributes.ID:
        id_info=entry.attributes.ID[0]
    else:
        print(entry.attributes)
        id_info=entry.attributes
    seq = [seq for seq in sequence_list if seq.id == entry.seqid]
    if test_for_multiple_gff_matches(seq, entry):
        new_seq = 'NA'
    else:
        seq_start = entry.start - 1
        seq_end = entry.end
        new_seq = SeqRecord(
            seq[0].seq[seq_start:seq_end],
            id=id_info,
            description=entry.attributes.Note[0]
            )
        if entry.strand == '-':
            rc_seq = new_seq.seq.reverse_complement()
            new_seq.seq = rc_seq
    return(new_seq)


def fix_id_info(entry):
    if type(entry.attributes) != str:
        id_info=entry.attributes.ID[0]
    else:
        #print(entry.attributes)
        id_info=entry.attributes
    return(id_info)


def extend_CDS_record(entry, id_info, seq, new_seq):
    exon_start = entry.start -1
    exon_end = entry.end
    new_exon = SeqRecord(
        seq[0].seq[exon_start:exon_end],
        id=id_info)
    if entry.strand == '-':
        rc_seq = new_exon.seq.reverse_complement()
        new_exon.seq = rc_seq
    extended_CDS = new_seq.seq + new_exon.seq
    final_seq = new_seq
    final_seq.seq = extended_CDS
    return(final_seq)
    


def fix_name(id_info, gene_gff_list):
    if '-RA' in id_info:
        gene = [gene for gene in gene_gff_list if id_info.split('-RA')[0] == gene.attributes.ID[0]]
    else:
        gene = [gene for gene in gene_gff_list if id_info == gene.attributes.ID[0]]
    return(gene)


def add_seq_description(entry, gene, exon_start, exon_end, id_info, seq):
    if type(entry.attributes) == 'list' and len(gene[0].attributes.Note) > 0 :
        new_seq = SeqRecord(
            seq[0].seq[exon_start:exon_end],
            id=id_info,
            description=gene[0].attributes.Note[0]
            )
    else:
        new_seq = SeqRecord(
            seq[0].seq[exon_start:exon_end],
            id=id_info
            )
    return(new_seq)


def initiate_CDS_record(entry, id_info, gene_gff_list, seq):
    gene = fix_name(id_info, gene_gff_list)
    if test_for_multiple_gff_matches(seq, entry) :
        new_seq = 'NA'
    else:
        #CDS_ID = id_info
        exon_start = entry.start -1
        exon_end = entry.end
        new_seq = add_seq_description(entry, gene, exon_start, exon_end, id_info, seq)
        if entry.strand == '-':
            rc_seq = new_seq.seq.reverse_complement()
            new_seq.seq = rc_seq
    return(new_seq)



def pull_cds(sequence_list, gff_list) :
    gene_gff_list = [ entry for entry in gff_list if entry.type == 'gene']
    CDS_gff_list = [ entry for entry in gff_list if entry.type == 'CDS']       
    feature_seqs = []
    CDS_ID = ''
    for entry in CDS_gff_list:
        id_info=fix_id_info(entry)
        if id_info != CDS_ID:
            seq = [seq for seq in sequence_list if seq.id == entry.seqid]
            new_seq = initiate_CDS_record(entry, id_info, gene_gff_list, seq)
            if not isinstance(new_seq, str) :  
                feature_seqs.append(new_seq)
            else:
                print('hey, an NA')  
        else:
            new_seq = extend_CDS_record(entry, id_info, seq, new_seq)
            if not isinstance(new_seq, str) :  
                feature_seqs[-1] = new_seq
            else:
                print('hey, an NA')  
        CDS_ID = id_info
    return(feature_seqs)


def Pull_gene_features(sequence_list, gff_list, feature) :
    if feature == 'gene':
        feature_seqs = pull_gene(sequence_list, gff_list, feature)
    elif feature == 'CDS':
        feature_seqs = pull_cds(sequence_list, gff_list)
    return(feature_seqs)


def remove_cds_string_in_id (seqlist):
    for seq in seqlist:
        if ':cds' in seq.id:
            seq.id = seq.id.split(':cds')[0]
    return(seqlist)


########################################
################# CODE #################
########################################

print("Reading gff")
gff_list = GFF_parse(gff)

print("Reading genome sequences")		
sequences = list(SeqIO.parse(seqfile,"fasta"))

print("Matching features with sequences")
feature_seqs = Pull_gene_features(sequences, gff_list, feature)
print("Recovered " + str(len(feature_seqs)) + " feature sequences")

if remove_cds_label == True :
    remove_cds_string_in_id (feature_seqs)
    

counter = 1
for seq in feature_seqs:
    if seq.id == 'unnamed':
        new_name = 'unnamed' + str(counter)
        seq.id = new_name
        counter += 1

print("Writing feature sequences to fasta")
SeqIO.write(feature_seqs, output, "fasta")

print("Finished.")
