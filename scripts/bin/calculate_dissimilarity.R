
library(argparser)
library(stringr)
library(vegan)
library(tibble)
library(hash)
library(ggplot2)
library(cowplot)
library(BiocGenerics)
library(gggenes)
library(grid)

print(1)

parser <- arg_parser("Rscript to extract dissimiliarity information")
# specify our desired options 
# by default ArgumentParser will add an help option 
parser<-add_argument(parser, "--orthogroup", type="character",
                    help="The focal orthogroup")
parser<-add_argument(parser, "--genomes_dir", type="character",default = "0_Inputs/genomes",
                     help="Directory containing genomes")
#parser<-add_argument(parser, "--species", type="character",
#                     help="Comma separated pair of species. Example: Species1,Species2")
#parser<-add_argument(parser, "--samples", type="character", 
#                    help="Comma separated pair of samples (one from each species) that can be found in the featureCounts output. Example: Sample1,Sample2. Can also use 'Average' for one or both species.")
parser<-add_argument(parser, "--aminoacid_dir", type="character", default = "3.7_OG_AA_kmers",
                    help="Directory where there are amino acid kmer counts for the orthogroup. default is 3.7_OG_AA_kmers")
parser<-add_argument(parser, "--orthogroups_info_dir", type="character", default = "3.6_OG_classifications",
                    help="Directory where orthogroup info and classifications are. default is 3.6_OG_classifications")
parser<-add_argument(parser, "--expression_dir", type="character", default = "3.8_OG_expression",
                    help="Directory where orthogroup info and classifications are. default is 3.8_OG_expression")
parser<-add_argument(parser, "--geneorder_dir", type="character", default = "3.9_OG_geneorder",
                     help="Directory containing a csv file of gene positions. default is 3.9_OG_geneorder")
parser<-add_argument(parser, "--output_dir", type="character", default = "4.1_OG_Routputs",
                    help="Directory where to write output files. default is 4.1_OG_Routputs")
args <- parse_args(parser)


print(2)

##species = "Scatenatus,Stergeminus"
##samples = "Scate-DRR0043-VG,Sterg-CLP2594-VG"

orthogroup = args$orthogroup
species = args$species
samples = args$samples
aminoacid_dir = args$aminoacid_dir
orthogroups_info_dir = args$orthogroups_info_dir
expression_dir = args$expression_dir
geneorder_dir = args$geneorder_dir
output_dir = args$output_dir


################################################################################################


## Count kmers
count_kmers = function(kmers_object){
  kmer_counts<-as.data.frame(table(unique(kmers_object$kmer)))
  return(kmer_counts)
}

print(c("Test1"))
## calculate kmer diversity
kmer_diversity = function(kmer_counts){
  ### transpose the data to calculate H
  print(ncol(kmer_counts))
  count_vector<-kmer_counts[,2]
  transposed_kmers<-t(count_vector)
  colnames(transposed_kmers) = kmer_counts[,1]
  
  ### calculate H as effective number of kmers
  H <- exp(diversity(transposed_kmers))
  #### Format output entry
  print("Guess it worked?")
  return(H)
}

print(c("Test2"))

## Summarize kmer diversity
Change_in_kmer_diversity = function(kmers_object){
  kmer_contribution_output =data.frame(Seq=character(),Total_H=numeric(),New_H=numeric(),
                                       H_change=numeric())
  ### create list of sequences
  seq_list<-unique(kmers_object$seqid)
  
  Total_H<-kmer_diversity(count_kmers(kmers_object))
  
  ## for each sequence, remove those kmers, recount H, calculate difference
  print(length(seq_list))
  if (length(seq_list) > 1){
    print("went down one? what?")
  for (seq in seq_list){
    other_kmers=subset(kmers_object,kmers_object$seqid != seq)
    counts<-count_kmers(other_kmers)
    New_H<-kmer_diversity(counts)
    print("Do we get here? and more than once?")
    H_change<-abs(Total_H-New_H)
    row<-tibble_row(seq,Total_H,New_H,H_change)
    kmer_contribution_output <- rbind(kmer_contribution_output,row)
  }} else{
    New_H<-0
    H_change<-abs(Total_H-New_H)
    row<-tibble_row(seq_list[1],Total_H,New_H,H_change)
    kmer_contribution_output <- rbind(kmer_contribution_output,row)
  } 
  
  return(kmer_contribution_output)
}



print(c("Test3"))

Calc_kmer_similarity<-function(kmers_data1,kmers_data2){
  ## Scate kmer table
  kmer_table1<-as.data.frame(table(kmers_data1$kmer))
  ## Sterg kmer table
  kmer_table2<-as.data.frame(table(kmers_data2$kmer))
  
  ## Scate_uniq_kmers
  table1_uniq_kmers<-subset(kmer_table1$Var1, !(kmer_table1$Var1 %in% kmer_table2$Var1))
  ## Sterg_uniq_kmers
  table2_uniq_kmers<-subset(kmer_table2$Var1, !(kmer_table2$Var1 %in% kmer_table1$Var1))
  
  ##Scate_additions
  table1_additions<-cbind.data.frame(table2_uniq_kmers,rep(0,length(table2_uniq_kmers)))
  colnames(table1_additions)<-c("Var1","Freq")
  ##Sterg_additions
  table2_additions<-cbind.data.frame(table1_uniq_kmers,rep(0,length(table1_uniq_kmers)))
  colnames(table2_additions)<-c("Var1","Freq")
  
  kmer_table1<-rbind(kmer_table1,table1_additions)
  colnames(kmer_table1)<-c("kmer","Species1")
  kmer_table2<-rbind(kmer_table2,table2_additions)
  colnames(kmer_table2)<-c("kmer","Species2")
  
  Species1_Species2_kmer_table<-merge(kmer_table1,kmer_table2, by.x="kmer")
  row.names(Species1_Species2_kmer_table)<-Species1_Species2_kmer_table$kmer
  Species1_Species2_kmer_table<-Species1_Species2_kmer_table[,2:3]
  Species1_Species2_kmer_table<-t(Species1_Species2_kmer_table)
  
  Species1_Species2_similarity<-vegdist(Species1_Species2_kmer_table)
  
  return(Species1_Species2_similarity)
}


print(c("Test4"))

##jacknife_orthogroups<-function(kmers_data1,kmers_data2,ortho_class_hash){
##  output<-data.frame(orthogroup=factor(),class=factor(),similarity=numeric())
##  all_orthogroups<-c(as.character(kmers_data1[,1]),as.character(kmers_data2[,1]))
##  #print(all_orthogroups)
##  uniq_orthos<-unique(all_orthogroups)
##  for (ortho in uniq_orthos){
##    reduced_kmers_1<-subset(kmers_data1,kmers_data1[,1] != ortho)
##    reduced_kmers_2<-subset(kmers_data2,kmers_data2[,1] != ortho)
##    similarity<-Calc_kmer_similarity(reduced_kmers_1,reduced_kmers_2)
##    class<-hash::values(ortho_class_hash,ortho)
##    entry<-tibble(ortho,class,similarity)
##    output<-rbind(output,entry)
##  }
##  return(output)
##}



Make_kmer_expression_table<-function(kmer_expr_object){
  kmer_tab<-data.frame(kmer=character(),expression=numeric())
  colnames(kmer_tab)<-c("kmer","expression")
  for (kmer_string in unique(kmer_expr_object$kmer)){
    #for (kmer_string in list){
    expression_kmer_entry<-subset(kmer_expr_object,kmer_expr_object$kmer == kmer_string)
    expression<-sum(expression_kmer_entry$expression)
    entry<-tibble(kmer_string, expression)
    colnames(entry)<-c("kmer","expression")
    kmer_tab<-rbind(kmer_tab,entry)
  } 
  return(kmer_tab)
}


print(c("Test5"))

Calc_kmer_expr_similarity<-function(kmer_expr_tab1,kmer_expr_tab2){
  ## Scate kmer table
  
  ## Scate_uniq_kmers
  table1_uniq_kmers<-subset(kmer_expr_tab1$kmer, !(kmer_expr_tab1$kmer %in% kmer_expr_tab2$kmer))
  
  ## Sterg_uniq_kmers
  table2_uniq_kmers<-subset(kmer_expr_tab2$kmer, !(kmer_expr_tab2$kmer %in% kmer_expr_tab1$kmer))
  
  ##Scate_additions
  table1_additions<-cbind.data.frame(table2_uniq_kmers,rep(0,length(table2_uniq_kmers)))
  colnames(table1_additions)<-c("kmer","expression")
  
  ##Sterg_additions
  table2_additions<-cbind.data.frame(table1_uniq_kmers,rep(0,length(table1_uniq_kmers)))
  colnames(table2_additions)<-c("kmer","expression")
  
  kmer_expr_tab1<-rbind(kmer_expr_tab1,table1_additions)
  colnames(kmer_expr_tab1)<-c("kmer","Species1_expression")
  kmer_expr_tab2<-rbind(kmer_expr_tab2,table2_additions)
  colnames(kmer_expr_tab2)<-c("kmer","Species2_expression")
  
  Species1_Species2_kmer_expr_table<-merge(kmer_expr_tab1,kmer_expr_tab2, by.x="kmer")
  row.names(Species1_Species2_kmer_expr_table)<-Species1_Species2_kmer_expr_table$kmer
  Species1_Species2_kmer_expr_table_expression<-Species1_Species2_kmer_expr_table[,c(2,3)]
  #Species1_Species2_kmer_expr_table_count<-Species1_Species2_kmer_expr_table[,c(2,4)]
  Species1_Species2_kmer_expr_table_expression<-t(Species1_Species2_kmer_expr_table_expression)
  
  Species1_Species2_similarity_expression<-vegdist(Species1_Species2_kmer_expr_table_expression)
  Species1_Species2_similarity<-c(Species1_Species2_similarity_expression)
  
  return(Species1_Species2_similarity)
}



define_seq_groups <- function(gene_order){
  group_list <-list()
  seqs = unique(gene_order$seq)
  for (i in 1:length(seqs)){
    group<-c()
    ref_subset <- subset(gene_order,gene_order$seq == seqs[i])
    for (j in 1:length(seqs)){
      test_subset <- subset(gene_order,gene_order$seq == seqs[j])
      if (any(test_subset$gene %in% ref_subset$gene)){
        group<-c(group,as.character(seqs[j]))}
      }  
    group_list[[i]]<-group
    }
  
  merged_groups<-list()
  merged_groups[1]<-group_list[1]
  for (j in 1:length(merged_groups)){
    for (i in 1:length(group_list)){
      if (any(group_list[[i]] %in% merged_groups[[j]])){
        merged_groups[j] <- unique(list(group_list[[i]],merged_groups[[j]]))
      } else{merged_groups[length(merged_groups)+1] <- group_list[i]}
    }  
  }
  seq_dict<-hash()
  counter = 0
  for (group in merged_groups){
    counter <- counter +1
    val<-paste("seq.",counter,sep = '')
    for (seq in group) {
      seq_dict[seq]<-val
    }
  }
  
  seq_groups<-hash::values(seq_dict,gene_order$seq)
  gene_order<-cbind(gene_order,seq_groups)
  return(gene_order)
}


print(c("Test6"))
################################################################################################
## Section for making plots of the genes/gene orders
gene_order_file = paste(geneorder_dir,"/",orthogroup,"/",orthogroup,"_positions.csv", sep='')
gene_order<-read.csv(gene_order_file)

#gene_order<-gene_order[-c(14),]

#dummies<-make_alignment_dummies(gene_order, aes(xmin = start, xmax = end, y = molecule, id = gene), on = gene_order$gene[1])

gene_order<-define_seq_groups(gene_order)


genes_plot<-ggplot(gene_order, aes(xmin = start, xmax = end, y = molecule, fill = gene, forward=orientation)) +
  geom_gene_arrow() +
  #geom_blank(data = dummies) +
  facet_grid(species ~ seq_groups, scales = "free", space = "free") +
  #scale_fill_brewer(palette = "Set3") +
  theme_genes()


plot_file = paste(output_dir,"/",orthogroup,"/",orthogroup,"_genes_plot.png",sep="")
print(genes_plot)
ggsave(plot_file,genes_plot, device="png")


################################################################################################
## define species strings
species1 = str_split(species, ",")[[1]][1]
species2 = str_split(species, ",")[[1]][2]

##print(c("Test7"))

## define sample strings
sample1 = str_split(samples, ",")[[1]][1]
sample2 = str_split(samples, ",")[[1]][2]
sample1 = str_replace_all(sample1, '-', '.')
sample2 = str_replace_all(sample2, '-', '.')


##print(c("Test8"))

## read in files
species_1_file = paste(aminoacid_dir,"/",orthogroup,"/",species1,"_",orthogroup,"_AA_kmers.csv",sep='')
species1_OG_AA_kmers <- read.csv(species_1_file)

species_2_file = paste(aminoacid_dir,"/",orthogroup,"/",species2,"_",orthogroup,"_AA_kmers.csv",sep='')
species2_OG_AA_kmers <- read.csv(species_2_file)

print(c("Test9"))
## calculate change in kmer diversity for each species
print(ncol(species1_OG_AA_kmers))
print(ncol(species2_OG_AA_kmers))
species1_change_kmer_div<-Change_in_kmer_diversity(species1_OG_AA_kmers)
print("So we pass Species 1")
species2_change_kmer_div<-Change_in_kmer_diversity(species2_OG_AA_kmers)
print("okay, pass Species 2")

print(1)
## Need to use our gene map and gene classes to create a classes list for each species
ortho_keys_file = paste(orthogroups_info_dir,"/",orthogroup,"/",orthogroup,"_genes.csv",sep='')
OG_ortho_key <- read.csv(ortho_keys_file)

print("K so we read the ortho_keys_file")

ortho_keys<-hash()
for (i in 1:nrow(OG_ortho_key)){
  ortho_keys[OG_ortho_key$gene[i]]<-OG_ortho_key$orthogroup[i]
}

print("and made the orthokeys hash")

OG_ortho_classes_file <- paste(orthogroups_info_dir,"/",orthogroup,"/",orthogroup,"_classes.csv",sep='')
OG_ortho_classes<-read.csv(OG_ortho_classes_file)

print("we read the orthoclasses too")

## select only the relevant two species colums
possible_string1 = paste(species1,species2,sep=".")
possible_string2 = paste(species2,species1,sep=".")
if (possible_string1 %in% colnames(OG_ortho_classes)){
  OG_classes<-as.data.frame(cbind(OG_ortho_classes$paralog,OG_ortho_classes[possible_string1]))
} else if (possible_string2 %in% colnames(OG_ortho_classes)){
  OG_classes<-as.data.frame(cbind(OG_ortho_classes$paralog,OG_ortho_classes[possible_string2]))
}
colnames(OG_classes)<-c("OG_orthos","OG_classes")
OG_classes<-na.omit(OG_classes)

print("made the classes table")

ortho_classes<-hash()
for (i in 1:length(OG_classes$OG_orthos)){
  ortho_classes[OG_classes$OG_orthos[i]]<-OG_classes$OG_classes[i]
}

print("made the classes hash")

species1_change_kmer_div<-as.data.frame(species1_change_kmer_div)
species1_change_kmer_div
print(colnames(species1_change_kmer_div))
print(is.data.frame(species1_change_kmer_div))
print(species1_change_kmer_div[["seq"]])
species1_gene_classes<-as.vector(hash::values(ortho_classes,as.vector(hash::values(ortho_keys,species1_change_kmer_div[["seq"]]))))

species2_change_kmer_div<-as.data.frame(species2_change_kmer_div)
print(colnames(species2_change_kmer_div))
species2_gene_classes<-as.vector(hash::values(ortho_classes,as.vector(hash::values(ortho_keys,species2_change_kmer_div$seq))))

print("made gene class vectors")
  

species1_change_kmer_div<-cbind(species1_change_kmer_div,species1_gene_classes)
species2_change_kmer_div<-cbind(species2_change_kmer_div,species2_gene_classes)


species1_change_perc<-species1_change_kmer_div$H_change/species1_change_kmer_div$Total_H*100
species1_change_kmer_div<-cbind(species1_change_kmer_div,species1_change_perc)
species2_change_perc<-species2_change_kmer_div$H_change/species2_change_kmer_div$Total_H*100
species2_change_kmer_div<-cbind(species2_change_kmer_div,species2_change_perc)

colnames(species1_change_kmer_div)<-c("seq","Total_H","New_H","H_change","gene_classes","change_perc")
colnames(species2_change_kmer_div)<-c("seq","Total_H","New_H","H_change","gene_classes","change_perc")

All_species_change_kmer_div<-rbind(species1_change_kmer_div,species2_change_kmer_div)

species_labels<-c(rep("species1",nrow(species1_change_kmer_div)),rep("species2",nrow(species2_change_kmer_div)))

All_species_change_kmer_div<-cbind(All_species_change_kmer_div,species_labels)

print(2)
## our seq_list_3 is OG_ortho_key$gene
AA_div_orthogroups<-hash::values(ortho_keys, All_species_change_kmer_div$seq)
All_species_change_kmer_div<-cbind(All_species_change_kmer_div,AA_div_orthogroups)

p<-ggplot(All_species_change_kmer_div, aes(x=AA_div_orthogroups,y=change_perc,fill=species_labels))+
geom_bar(stat="identity", position = position_dodge()) + theme_cowplot() +
  theme(axis.text.x = element_text(angle = 90)) + scale_fill_discrete(labels=c(species1, species2))
p



## Modify the Scate and Sterg Kmer tables to have orthogroup names
#species1_kmer_seqs<-c()
#for (i in 1:nrow(species1_OG_AA_kmers)){
#  SVMP<-str_replace(as.character(species1_OG_AA_kmers$seqid[i]),':cds','')
#  species1_kmer_seqs<-c(species1_kmer_seqs,species1_OG_AA_kmers$seqid)
#}
species1_ortho_assignments<-hash::values(ortho_keys,species1_OG_AA_kmers$seqid)


#Sterg_kmer_seqs<-c()
#for (i in 1:nrow(Stergeminus_SVMP_AA_kmers)){
#  SVMP<-str_replace(as.character(Stergeminus_SVMP_AA_kmers$seqid[i]),':cds','')
#  Sterg_kmer_seqs<-c(Sterg_kmer_seqs,SVMP)
#}
species2_ortho_assignments<-hash::values(ortho_keys,species2_OG_AA_kmers$seqid)


OG_cats<-hash()
for (i in 1:nrow(All_species_change_kmer_div)){
  OG_cats[All_species_change_kmer_div$seq[i]]<-All_species_change_kmer_div$gene_classes[i]
}

species1_class_assignments<-hash::values(OG_cats,species1_OG_AA_kmers$seqid)

species2_class_assignments<-hash::values(OG_cats,species2_OG_AA_kmers$seqid)


species1_reorg_kmers<-cbind(species1_ortho_assignments,species1_class_assignments,
                         species1_OG_AA_kmers)

species2_reorg_kmers<-cbind(species2_ortho_assignments,species2_class_assignments,
                         species2_OG_AA_kmers)



print(3)
##kmer_jacknife<-jacknife_orthogroups(Scate_reorg_Kmers,Sterg_reorg_Kmers,ortho_classes)
##kmer_jacknife$similarity<-as.numeric(as.character(kmer_jacknife$similarity))

##ct_kmer_orthologs<-subset(kmer_jacknife,kmer_jacknife$class == "conserved_ortholog" | kmer_jacknife$class == "ct-paralog")
##ortholog_list<-rep("orthologs",nrow(ct_kmer_orthologs))
##ct_kmer_orthologs<-cbind(ct_kmer_orthologs,ortholog_list)
##colnames(ct_kmer_orthologs)<-c("ortho","class","dissimilarity","ortho_para")

##ct_kmer_paralogs<-subset(kmer_jacknife,kmer_jacknife$class != "conserved_ortholog" & kmer_jacknife$class != "ct-paralog")
##paralog_list<-rep("paralogs",nrow(ct_kmer_paralogs))
##ct_kmer_paralogs<-cbind(ct_kmer_paralogs,paralog_list)
##colnames(ct_kmer_paralogs)<-c("ortho","class","dissimilarity","ortho_para")


##ct_kmer_table<-rbind(ct_kmer_orthologs,ct_kmer_paralogs)

#Scate_Sterg_similarity<-Calc_kmer_similarity(Scate_reorg_Kmers,Sterg_reorg_Kmers)

#ggboxplot(ct_kmer_table, x="ortho_para", y="dissimilarity",add = "dotplot", fill = "cyan3", add.params=list(size=0.6),xlab = "Removed toxin type")+font("x.text", size = 10)+ geom_hline(yintercept=Scate_Sterg_similarity, linetype="dashed", color = "red", size=2)

##ggboxplot(ct_kmer_paralogs, x="class", y="dissimilarity",add = "dotplot", fill = "cyan3", add.params=list(size=0.6),xlab = "Removed toxin type")+font("x.text", size = 10)+ geom_hline(yintercept=Scate_Sterg_similarity, linetype="dashed", color = "red", size=2)


## Overall similarity
species1_species2_similarity<-Calc_kmer_similarity(species1_reorg_kmers,species2_reorg_kmers)

## Minus conserved orthologs
species1_paralog_kmers<-subset(species1_reorg_kmers,species1_reorg_kmers$species1_class_assignments != "conserved")
species2_paralog_kmers<-subset(species2_reorg_kmers,species2_reorg_kmers$species2_class_assignments != "conserved")
if (nrow(species1_paralog_kmers) == 0 | nrow(species1_paralog_kmers) == 0){paralog_similarity<-NA} else{
  paralog_similarity<-Calc_kmer_similarity(species1_paralog_kmers,species2_paralog_kmers)  
}


## Minus paralogs (all kinds)
species1_orthos_kmers<-subset(species1_reorg_kmers,species1_reorg_kmers$species1_class_assignments == "conserved")
species2_orthos_kmers<-subset(species2_reorg_kmers,species2_reorg_kmers$species2_class_assignments == "conserved" )
if (nrow(species1_orthos_kmers) == 0 | nrow(species1_orthos_kmers) == 0){orthos_similarity<-NA} else{
  orthos_similarity<-Calc_kmer_similarity(species1_orthos_kmers,species2_orthos_kmers)
}

## Minus in-paralogs
#Scate_noninparalog_kmers<-subset(Scate_reorg_Kmers,Scate_reorg_Kmers$Scate_class_assignments != "in-paralog")
#Sterg_noninparalog_kmers<-subset(Sterg_reorg_Kmers,Sterg_reorg_Kmers$Sterg_class_assignments != "in-paralog")
#in_paralog_similarity<-Calc_kmer_similarity(Scate_noninparalog_kmers,Sterg_noninparalog_kmers)


## Minus lost-paralogs
#Scate_nonlostparalog_kmers<-subset(Scate_reorg_Kmers,Scate_reorg_Kmers$Scate_class_assignments != "ca-paralog")
#Sterg_nonlostparalog_kmers<-subset(Sterg_reorg_Kmers,Sterg_reorg_Kmers$Sterg_class_assignments != "ca-paralog")
#lost_paralog_similarity<-Calc_kmer_similarity(Scate_nonlostparalog_kmers,Sterg_nonlostparalog_kmers)


dissimilarities<-c(species1_species2_similarity[1],paralog_similarity[1],
                   orthos_similarity[1])
#dissimilarities<-c(species1_species2_similarity[1],paralog_similarity[1],
#                   orthos_similarity[1],in_paralog_similarity[1], 
#                   lost_paralog_similarity[1])
groups<-c("Overall","Paralogs","Orthologs")
#groups<-c("Overall","Paralogs","Orthologs","In-paralogs","Out-paralogs")
expr_class<-rep("kmers_only",length(dissimilarities))

kmer_group_similarities<-cbind.data.frame(groups,dissimilarities,expr_class)

ortho_all_similarity<-subset(kmer_group_similarities, kmer_group_similarities$groups == "Orthologs" | kmer_group_similarities$groups == "Overall")

p<-ggplot(data=ortho_all_similarity, aes(x=groups, y=dissimilarities)) +
  geom_bar(stat="identity", fill=c("cyan3")) + labs(x="Gene Class", y="Dissimilarity")
p

print(4)
#paralog_similarity_tab<-subset(kmer_group_similarities, kmer_group_similarities$groups == "In-paralogs" | kmer_group_similarities$groups == "Out-paralogs" | kmer_group_similarities$groups == "Overall")

#p<-ggplot(data=paralog_similarity_tab, aes(x=groups, y=dissimilarities)) +
#  geom_bar(stat="identity", fill=c("green","red", "blue")) + labs(x="Gene Class", y=" Disimilarity")
#p



## First, just add expression counts to the kmer table
## First, step create hashes
OG_expression_file  = paste(expression_dir,"/",orthogroup,"/",orthogroup,"_expression.csv",sep="")
OG_expression = read.csv(OG_expression_file)

print(5)

if (sample1 %in% colnames(OG_expression) == FALSE){print("Uh-oh. Looks like your sample isnt in your expression file")}

if (sample2 %in% colnames(OG_expression) == FALSE){print("Uh-oh. Looks like your sample isnt in your expression file")}


species1_expression_hash<-hash()
for (i in 1:nrow(OG_expression)){
  species1_expression_hash[OG_expression$Orthogroup[i]]<-OG_expression[sample1][[1]][i]
}

species2_expression_hash<-hash()
for (i in 1:nrow(OG_expression)){
  species2_expression_hash[OG_expression$Orthogroup[i]]<-OG_expression[sample2][[1]][i]
}

print(6)

## Then generate the expression vectors
species1_expression_count<-hash::values(species1_expression_hash,species1_reorg_kmers$species1_ortho_assignments)
species2_expression_count<-hash::values(species2_expression_hash,species2_reorg_kmers$species2_ortho_assignments)


species1_kmer_expr<-cbind(species1_reorg_kmers,species1_expression_count)
colnames(species1_kmer_expr)<-c("species1_ortho_assignments","species1_class_assignments","seqid",
                             "kmer","expression")
species2_kmer_expr<-cbind(species2_reorg_kmers,species2_expression_count)
colnames(species2_kmer_expr)<-c("species2_ortho_assignments","species2_class_assignments","seqid",
                             "kmer","expression")



species1_kmer_expr_tab<-Make_kmer_expression_table(species1_kmer_expr)
species1_kmer_expr_tab<-na.omit(species1_kmer_expr_tab)

species2_kmer_expr_tab<-Make_kmer_expression_table(species2_kmer_expr)
species2_kmer_expr_tab<-na.omit(species2_kmer_expr_tab)



species1_species2_expr_similarity<-Calc_kmer_expr_similarity(species1_kmer_expr_tab,species2_kmer_expr_tab)[1]


print(7)

#jacknife_orthogroups_expr<-function(kmers_expr1,kmers_expr2,ortho_class_hash){
#  output<-data.frame(orthogroup=factor(),class=factor(),similarity=numeric())
#  all_orthogroups<-c(as.character(kmers_expr1[,1]),as.character(kmers_expr2[,1]))
#  #print(all_orthogroups)
#  uniq_orthos<-unique(all_orthogroups)
#  for (ortho in uniq_orthos){
#    red_kmers_expr_1<-subset(kmers_expr1,kmers_expr1[,1] != ortho)
#    red_kmers_expr_2<-subset(kmers_expr2,kmers_expr2[,1] != ortho)
#    kmers_expr_tab1<-Make_kmer_expression_table(red_kmers_expr_1)
#    kmers_expr_tab2<-Make_kmer_expression_table(red_kmers_expr_2)
#    similarity<-Calc_kmer_expr_similarity(kmers_expr_tab1,kmers_expr_tab2)[1]
#    class<-hash::values(ortho_class_hash,ortho)
#    entry<-tibble(ortho,class,similarity)
#    output<-rbind(output,entry)
#  }
#  return(output)
#}


#kmer_expr_jacknife<-jacknife_orthogroups_expr(species1_kmer_expr,species2_kmer_expr,ortho_classes)




#species1_species2_kmer_expr_orthologs<-subset(kmer_expr_jacknife,kmer_expr_jacknife$class == "conserved")
#ortholog_list<-rep("orthologs",nrow(species1_species2_kmer_expr_orthologs))
#species1_species2_kmer_expr_orthologs<-cbind(species1_species2_kmer_expr_orthologs,ortholog_list)
#colnames(species1_species2_kmer_expr_orthologs)<-c("ortho","class","dissimilarity","ortho_para")

#species1_species2_kmer_expr_paralogs<-subset(kmer_expr_jacknife,kmer_expr_jacknife$class != "conserved")
#paralog_list<-rep("paralogs",nrow(species1_species2_kmer_expr_paralogs))
#species1_species2_kmer_expr_paralogs<-cbind(species1_species2_kmer_expr_paralogs,paralog_list)
#colnames(species1_species2_kmer_expr_paralogs)<-c("ortho","class","dissimilarity","ortho_para")

#species1_species2_kmer_expr_table<-rbind(ct_kmer_expr_orthologs,ct_kmer_expr_paralogs)

#ggboxplot(species1_species2_kmer_expr_table, x="ortho_para", y="dissimilarity",add = "dotplot", fill = "indianred1", add.params=list(size=0.6),xlab = "Removed toxin type")+font("x.text", size = 10)+ geom_hline(yintercept=Scate_Sterg_expr_similarity, linetype="dashed", color = "red", size=2)

#ggboxplot(species1_species2_kmer_expr_paralogs, x="class", y="dissimilarity",add = "dotplot", fill = "indianred1", add.params=list(size=0.6),xlab = "Removed toxin type")+font("x.text", size = 10)+ geom_hline(yintercept=Scate_Sterg_expr_similarity, linetype="dashed",color = "red", size=2)


#kmer_list<-rep("kmer",nrow(species1_species2_kmer_table))
#expr_list<-rep("expr",nrow(species1_species2_kmer_table))
#kmer_expr<-c(kmer_list,expr_list)
#kmer_expr_tab<-rbind(ct_kmer_table,ct_kmer_expr_table)
#kmer_expr_tab<-cbind(kmer_expr_tab,kmer_expr)

#ggboxplot(kmer_expr_tab, x="ortho_para", y="dissimilarity", fill="kmer_expr", add.params=list(size=0.6),xlab = "Removed toxin type")+font("x.text", size = 10)


#p


print(8)

## Minus conserved orthologs
species1_paralog_kmers_expr<-subset(species1_kmer_expr,species1_kmer_expr$species1_class_assignments != "conserved")
species2_paralog_kmers_expr<-subset(species2_kmer_expr,species2_kmer_expr$species2_class_assignments != "conserved")
species1_paralog_kmers_expr_tab<-Make_kmer_expression_table(species1_paralog_kmers_expr)
species1_paralog_kmers_expr_tab<-na.omit(species1_paralog_kmers_expr_tab)
species2_paralog_kmers_expr_tab<-Make_kmer_expression_table(species2_paralog_kmers_expr)
species2_paralog_kmers_expr_tab<-na.omit(species2_paralog_kmers_expr_tab)
paralog_similarity_expr<-Calc_kmer_expr_similarity(species1_paralog_kmers_expr_tab, species2_paralog_kmers_expr_tab)[1]

print(9)
## Minus paralogs (all kinds)
species1_orthos_kmers_expr<-subset(species1_kmer_expr,species1_kmer_expr$species1_class_assignments == "conserved")
species2_orthos_kmers_expr<-subset(species2_kmer_expr,species2_kmer_expr$species2_class_assignments == "conserved")
species1_orthos_kmers_expr_tab<-Make_kmer_expression_table(species1_orthos_kmers_expr)
species1_orthos_kmers_expr_tab<-na.omit(species1_orthos_kmers_expr_tab)
species2_orthos_kmers_expr_tab<-Make_kmer_expression_table(species2_orthos_kmers_expr)
species2_orthos_kmers_expr_tab<-na.omit(species2_orthos_kmers_expr_tab)
orthos_similarity_expr<-Calc_kmer_expr_similarity(species1_orthos_kmers_expr_tab,
                                                  species2_orthos_kmers_expr_tab)[1]

## Minus in-paralogs
##species1_noninparalog_kmers_expr<-subset(species1_kmer_expr,species2_kmer_expr$species1_class_assignments != "duplcation")
##species2_noninparalog_kmers_expr<-subset(species2_kmer_expr,species2_kmer_expr$species2_class_assignments != "duplication")
##species1_noninparalog_kmers_expr_tab<-Make_kmer_expression_table(species1_noninparalog_kmers_expr)
##species2_noninparalog_kmers_expr_tab<-Make_kmer_expression_table(species2_noninparalog_kmers_expr)
##in_paralog_similarity_expr<-Calc_kmer_expr_similarity(species1_noninparalog_kmers_expr_tab, species2_noninparalog_kmers_expr_tab)[1]


## Minus lost-paralogs
##species1_nonlostparalog_kmers_expr<-subset(species1_kmer_expr,species1_kmer_expr$species1_class_assignments != "loss")
##species2_nonlostparalog_kmers_expr<-subset(species2_kmer_expr,species2_kmer_expr$species2_class_assignments != "loss")
##species1_nonlostparalog_kmers_expr_tab<-Make_kmer_expression_table(species1_nonlostparalog_kmers_expr)
##species2_nonlostparalog_kmers_expr_tab<-Make_kmer_expression_table(species2_nonlostparalog_kmers_expr)
##lost_paralog_similarity_expr<-Calc_kmer_expr_similarity(species1_nonlostparalog_kmers_expr_tab,species2_nonlostparalog_kmers_expr_tab)[1]

print(10)

expr_dissimilarities<-c(species1_species2_expr_similarity[1],orthos_similarity_expr[1],
                        paralog_similarity_expr[1])
groups<-c("Overall","Orthologs","Paralogs")
expr_class<-rep("expressed",rep(length(groups)))

expr_group_similarities<-cbind.data.frame(groups,expr_dissimilarities,expr_class)
colnames(expr_group_similarities)<-c("groups","dissimilarities","expr_class")

all_group_similarities<-rbind(kmer_group_similarities,expr_group_similarities)


## Here is the orthos and overall plot
ortho_para_expr_similarity<-subset(all_group_similarities, all_group_similarities$groups == "Orthologs" | all_group_similarities$groups == "Overall")

p<-ggplot(data=ortho_para_expr_similarity, aes(x=groups, y=dissimilarities, fill=expr_class)) +
  geom_bar(stat="identity", position=position_dodge()) + labs(x="Gene Class", y="Dissimilarity")


make_directory_command = paste("mkdir ",output_dir,"/",orthogroup,sep="")
system(make_directory_command)


plot_file = paste(output_dir,"/",orthogroup,"/",orthogroup,"_dissimilarity_plot.png",sep="")
print(plot_file)
ggsave(plot_file,p, device="png")

classifications<-c()
for (i in 1:nrow(ortho_para_expr_similarity)){
  classification = paste(ortho_para_expr_similarity[i,1],ortho_para_expr_similarity[i,3],sep="_")
  classifications<-c(classifications, classification)
}

output_line = t(ortho_para_expr_similarity[,2])
colnames(output_line) = classifications


output_table_file = paste(output_dir,"/",orthogroup,"/",orthogroup,"_dissimilarity_table.csv",sep = '')
write.csv(output_line,output_table_file, row.names = FALSE)












