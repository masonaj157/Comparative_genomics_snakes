
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
library(svglite)

print(1)

parser <- arg_parser("Rscript to extract dissimiliarity information")
# specify our desired options 
# by default ArgumentParser will add an help option 
parser<-add_argument(parser, "--orthogroup", type="character",
                    help="The focal orthogroup")
parser<-add_argument(parser, "--genomes_dir", type="character",default = "0_Inputs/genomes",
                     help="Directory containing genomes")
parser<-add_argument(parser, "--samples", type="character", default = "0_Inputs/sample_assignment.csv",
                    help="File containing species assignment for each sample used in RNAseq mapping. Default is 0_inputs/sample_assignment.csv")
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
samples = args$samples
aminoacid_dir = args$aminoacid_dir
orthogroups_info_dir = args$orthogroups_info_dir
expression_dir = args$expression_dir
geneorder_dir = args$geneorder_dir
genomes_dir = args$genomes_dir
output_dir = args$output_dir


make_directory_command = paste("mkdir ",output_dir,"/",orthogroup,sep="")
system(make_directory_command)

################################################################################################


## Count kmers
count_kmers = function(kmers_object){
  kmer_counts<-as.data.frame(table(unique(kmers_object$kmer)))
  return(kmer_counts)
}

##print(c("Test1"))
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

##print(c("Test2"))

## Summarize kmer diversity
Change_in_kmer_diversity = function(kmers_object){
  kmer_contribution_output =data.frame(Seq=character(),Total_H=numeric(),New_H=numeric(),
                                       H_change=numeric())
  ### create list of sequences
  seq_list<-unique(kmers_object$seqid)
  
  Total_H<-kmer_diversity(count_kmers(kmers_object))
  
  ## for each sequence, remove those kmers, recount H, calculate difference
  #print(length(seq_list))
  if (length(seq_list) > 1){
    print("went down one? what?")
  for (seq in seq_list){
    other_kmers=subset(kmers_object,kmers_object$seqid != seq)
    counts<-count_kmers(other_kmers)
    New_H<-kmer_diversity(counts)
    #print("Do we get here? and more than once?")
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



#print(c("Test3"))

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


Calc_kmer_disimilarity_effect<-function(kmers_data1,kmers_data2){
  ## Scate kmer table
  kmer_table1<-as.data.frame(table(kmers_data1$kmer))
  ## Sterg kmer table
  kmer_table2<-as.data.frame(table(kmers_data2$kmer))
  
  ## Species1_uniq_kmers
  table1_uniq_kmers<-subset(kmer_table1$Var1, !(kmer_table1$Var1 %in% kmer_table2$Var1))
  ## Species2_uniq_kmers
  table2_uniq_kmers<-subset(kmer_table2$Var1, !(kmer_table2$Var1 %in% kmer_table1$Var1))
  
  ##Species1_additions
  table1_additions<-cbind.data.frame(table2_uniq_kmers,rep(0,length(table2_uniq_kmers)))
  colnames(table1_additions)<-c("Var1","Freq")
  ##Species2_additions
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
  
  Species1_Species2_disimilarity<-vegdist(Species1_Species2_kmer_table)
  kmer_differences<-abs(Species1_Species2_kmer_table[1,]-Species1_Species2_kmer_table[2,])
  kmer_totals<-Species1_Species2_kmer_table[1,]+Species1_Species2_kmer_table[2,]
  sum_kmer_diff<-sum(kmer_differences)
  sum_kmer_totals<-sum(kmer_totals)
  
  orthos<-unique(c(kmers_data1$species_ortho_assignments,kmers_data2$species_ortho_assignments))
  
  ## For each unique orthogroup we will subset to their kmers, then use our formula to calculate dissimilarity effect
  ## Total differences and Total numbers will be pulled from the original transposed kmers table.
  disimilarity_effects<-c()
  for (ortho in orthos){
    ortho_kmers_1<-kmers_data1[kmers_data1$species_ortho_assignments == ortho,]
    ortho_kmers_2<-kmers_data2[kmers_data2$species_ortho_assignments == ortho,]
    ortho_kmer_table1<-as.data.frame(table(ortho_kmers_1$kmer))
    ortho_kmer_table2<-as.data.frame(table(ortho_kmers_2$kmer))
    #print(ortho)
    #print(1)
    ## Species1_uniq_kmers
    ortho_table1_uniq_kmers<-subset(ortho_kmer_table1$Var1, !(ortho_kmer_table1$Var1 %in% ortho_kmer_table2$Var1))
    ## Species2_uniq_kmers
    ortho_table2_uniq_kmers<-subset(ortho_kmer_table2$Var1, !(ortho_kmer_table2$Var1 %in% ortho_kmer_table1$Var1))
    #print(2)
    ##Species1_additions
    ortho_table1_additions<-cbind.data.frame(ortho_table2_uniq_kmers,rep(0,length(ortho_table2_uniq_kmers)))
    if (nrow(ortho_table1_additions) != 0){
    colnames(ortho_table1_additions)<-c("Var1","Freq")
    ortho_kmer_table1<-rbind(ortho_kmer_table1,ortho_table1_additions)
    }
    ##Species2_additions
    ortho_table2_additions<-cbind.data.frame(ortho_table1_uniq_kmers,rep(0,length(ortho_table1_uniq_kmers)))
    if (nrow(ortho_table2_additions) != 0){
    colnames(ortho_table2_additions)<-c("Var1","Freq")
    ortho_kmer_table2<-rbind(ortho_kmer_table2,ortho_table2_additions)
    }
    #print(3)
    colnames(ortho_kmer_table1)<-c("kmer","Species1")
    colnames(ortho_kmer_table2)<-c("kmer","Species2")
    #print(4)
    ortho_Species1_Species2_kmer_table<-merge(ortho_kmer_table1,ortho_kmer_table2, by.x="kmer")
    row.names(ortho_Species1_Species2_kmer_table)<-ortho_Species1_Species2_kmer_table$kmer
    ortho_Species1_Species2_kmer_table<-ortho_Species1_Species2_kmer_table[,2:3]
    ortho_Species1_Species2_kmer_table<-t(ortho_Species1_Species2_kmer_table)
    #print(5)
    ortho_kmer_differences<-abs(ortho_Species1_Species2_kmer_table[1,]-ortho_Species1_Species2_kmer_table[2,])
    ortho_kmer_totals<-ortho_Species1_Species2_kmer_table[1,]+ortho_Species1_Species2_kmer_table[2,]
    ortho_sum_kmer_diff<-sum(ortho_kmer_differences)
    ortho_sum_kmer_totals<-sum(ortho_kmer_totals)
    #print(6)
    Numerator<-((1-(ortho_sum_kmer_diff/sum_kmer_diff)) * Species1_Species2_disimilarity[[1]])
    Denominator<-(1-(ortho_sum_kmer_totals/sum_kmer_totals))
    ortho_dissimilarity<-Numerator/Denominator
    Dissimilarity_effect<- Species1_Species2_disimilarity[[1]] - ortho_dissimilarity
    disimilarity_effects<-c(disimilarity_effects,Dissimilarity_effect)
  }
  
  ortho_dissimilarity_effects=data.frame(orthos,disimilarity_effects)
  return(ortho_dissimilarity_effects)
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
  kmer_tab <- subset(kmer_tab,kmer_tab$expression != 0)
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



Calc_kmer_expr_disimilarity_effects<-function(species1_kmer_expr, species2_kmer_expr, kmer_expr_tab1,kmer_expr_tab2){
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
  
  kmer_differences<-abs(Species1_Species2_kmer_expr_table_expression[1,]-Species1_Species2_kmer_expr_table_expression[2,])
  kmer_totals<-Species1_Species2_kmer_expr_table_expression[1,]+Species1_Species2_kmer_expr_table_expression[2,]
  sum_kmer_diff<-sum(kmer_differences)
  sum_kmer_totals<-sum(kmer_totals)
  
  orthos<-unique(c(species1_kmer_expr$species_ortho_assignments,species2_kmer_expr$species_ortho_assignments))
  
  ## For each unique orthogroup we will subset to their kmers, then use our formula to calculate dissimilarity effect
  ## Total differences and Total numbers will be pulled from the original transposed kmers table.
  disimilarity_effects<-c()
  for (ortho in orthos){
    #print(ortho)
    ortho_kmers_1<-species1_kmer_expr[species1_kmer_expr$species_ortho_assignments == ortho,]
    ortho_kmers_2<-species2_kmer_expr[species2_kmer_expr$species_ortho_assignments == ortho,]
    ortho_species1_kmer_expr_tab<-Make_kmer_expression_table(ortho_kmers_1)
    ortho_species1_kmer_expr_tab<-na.omit(ortho_species1_kmer_expr_tab)
    ortho_species2_kmer_expr_tab<-Make_kmer_expression_table(ortho_kmers_2)
    ortho_species2_kmer_expr_tab<-na.omit(ortho_species2_kmer_expr_tab)
    
    ## Species1_uniq_kmers
    ortho_table1_uniq_kmers<-subset(ortho_species1_kmer_expr_tab$kmer, !(ortho_species1_kmer_expr_tab$kmer %in% ortho_species2_kmer_expr_tab$kmer))
    ## Species2_uniq_kmers
    ortho_table2_uniq_kmers<-subset(ortho_species2_kmer_expr_tab$kmer, !(ortho_species2_kmer_expr_tab$kmer %in% ortho_species1_kmer_expr_tab$kmer))
    #print(2)
    ##Species1_additions
    ortho_table1_additions<-cbind.data.frame(ortho_table2_uniq_kmers,rep(0,length(ortho_table2_uniq_kmers)))
    if (nrow(ortho_table1_additions) != 0){
      colnames(ortho_table1_additions)<-c("kmer","expression")
      ortho_species1_kmer_expr_tab<-rbind(ortho_species1_kmer_expr_tab,ortho_table1_additions)
    }
    ##Species2_additions
    ortho_table2_additions<-cbind.data.frame(ortho_table1_uniq_kmers,rep(0,length(ortho_table1_uniq_kmers)))
    if (nrow(ortho_table2_additions) != 0){
      colnames(ortho_table2_additions)<-c("kmer","expression")
      ortho_species2_kmer_expr_tab<-rbind(ortho_species2_kmer_expr_tab,ortho_table2_additions)
    }
    #print(3)
    colnames(ortho_species1_kmer_expr_tab)<-c("kmer","Species1")
    colnames(ortho_species2_kmer_expr_tab)<-c("kmer","Species2")
    #print(4)
    ortho_Species1_Species2_kmer_table<-merge(ortho_species1_kmer_expr_tab,ortho_species2_kmer_expr_tab, by.x="kmer")
    row.names(ortho_Species1_Species2_kmer_table)<-ortho_Species1_Species2_kmer_table$kmer
    ortho_Species1_Species2_kmer_table<-ortho_Species1_Species2_kmer_table[,2:3]
    ortho_Species1_Species2_kmer_table<-t(ortho_Species1_Species2_kmer_table)
    #print(5)
    ortho_kmer_differences<-abs(ortho_Species1_Species2_kmer_table[1,]-ortho_Species1_Species2_kmer_table[2,])
    ortho_kmer_totals<-ortho_Species1_Species2_kmer_table[1,]+ortho_Species1_Species2_kmer_table[2,]
    ortho_sum_kmer_diff<-sum(ortho_kmer_differences)
    ortho_sum_kmer_totals<-sum(ortho_kmer_totals)
    #print(6)
    Numerator<-((1-(ortho_sum_kmer_diff/sum_kmer_diff)) * Species1_Species2_similarity_expression[[1]])
    Denominator<-(1-(ortho_sum_kmer_totals/sum_kmer_totals))
    ortho_dissimilarity<-Numerator/Denominator
    Dissimilarity_effect<- Species1_Species2_similarity_expression[[1]] - ortho_dissimilarity
    disimilarity_effects<-c(disimilarity_effects,Dissimilarity_effect)
  }
  
  ortho_dissimilarity_effects=data.frame(orthos,disimilarity_effects)
  return(ortho_dissimilarity_effects)

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
        merged_groups[[j]] <- unique(c(group_list[[i]],merged_groups[[j]]))
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


#print(c("Test6"))
################################################################################################
## Section for making plots of the genes/gene orders
##print(c("Test1"))
#gene_order_file = paste(geneorder_dir,"/",orthogroup,"/",orthogroup,"_positions.csv", sep='')
#gene_order<-read.csv(gene_order_file)


#gene_order_dat<-define_seq_groups(gene_order)
######gene_order_dat<-subset(gene_order_dat,gene_order_dat$seq_groups=="seq.1")

######dummies<-make_alignment_dummies(gene_order_dat, aes(xmin = start, xmax = end, y = molecule, id = gene), on = gene_order_dat$gene[1])


#genes_plot<-ggplot(gene_order_dat, aes(xmin = start, xmax = end, y = molecule, fill = gene, forward=gene_order_dat$orientation)) +
  #geom_gene_arrow() + #theme(legend.position = "none") +
  #####geom_blank(data = dummies) +
  #facet_grid(species ~ seq_groups, scales = "free", space = "free") +
  #scale_fill_brewer(palette = "Set3") +
  #theme_genes()


##plot_file = paste(output_dir,"/",orthogroup,"/",orthogroup,"_genes_plot.svg",sep="")
##print(genes_plot)
##ggsave(plot_file,genes_plot, device="svg")
##plot_file = paste(output_dir,"/",orthogroup,"/",orthogroup,"_genes_plot.png",sep="")
##print(genes_plot)
##ggsave(plot_file,genes_plot, device="png")


#print(c("Test2"))

################################################################################################
#print(c("Test3"))

##get list of species from genomes dir
species_files<-list.files(genomes_dir)
species<-c()
for (i in 1:length(species_files)){
  sp = str_split(species_files[i], ".fasta")[[1]][1]
  species<-c(species,sp)
}


## Read in expression file, tpm10k. This will be used to 1) Remove species without representation, 2) later remove transcripts without representation
OG_expression_file  = paste(expression_dir,"/",orthogroup,"/",orthogroup,"_tpm10k_expression.csv",sep="")
OG_expression = read.csv(OG_expression_file)
sample_assignment<-read.csv(samples)
##print(sample_assignment)
colnames(sample_assignment)

## First remove species
absent_species<-c()
for (sp in species){
  print(sp)
  tmp_samples <- subset(sample_assignment,sample_assignment$Species == sp)
  for (i in 1:nrow(tmp_samples)){
    sample = str_replace_all(tmp_samples[i,1], '-', '.')
    tmp_samples[i,1]<-sample
  }
  tmp_expression <- OG_expression[,tmp_samples$Sample]
  if (sum(tmp_expression) == 0){
    absent_species<-c(absent_species, sp)  
  }
}

species<-species[!(species %in% absent_species)]
## Species without any expression have now been removed and will not be represented in the species comparisons


## Makes species-specific lists of unexpressed k-mers
for (sp in species){
  tmp_samples <- subset(sample_assignment,sample_assignment$Species == sp)
  for (i in 1:nrow(tmp_samples)){
    sample = str_replace_all(tmp_samples[i,1], '-', '.')
    tmp_samples[i,1]<-sample
  }
  tmp_expression <- OG_expression[,tmp_samples$Sample]
  if(is.data.frame(tmp_expression)){
    row_means<-rowMeans(tmp_expression)
    mean_expression<-cbind(OG_expression[,1],row_means) 
  } else {mean_expression<-cbind(OG_expression[,1],tmp_expression) }
  unexpressed<-subset(mean_expression,mean_expression[,2] == 0)
  unexpressed<-c(unexpressed[,1])
  unexpressed_name<-paste(sp,"_unexpressed_transcripts",sep='')
  assign(unexpressed_name,unexpressed)
}


species_pairs<-data.frame(species1=character(),species2=character())
for (i in 1:(length(species)-1)){
  #print(i)
  for (j in (i+1):length(species)){
    #print(j)
    species_pair<-c(species[i],species[j])
    species_pairs<-rbind(species_pairs,species_pair)
  }
}
colnames(species_pairs)<-c("species1","species2")



#print(c("Test7"))


## Need to use our gene map and gene classes to create a classes list for each species
ortho_keys_file = paste(orthogroups_info_dir,"/",orthogroup,"/",orthogroup,"_genes.csv",sep='')
OG_ortho_key <- read.csv(ortho_keys_file)

#print("K so we read the ortho_keys_file")

ortho_keys<-hash()
for (i in 1:nrow(OG_ortho_key)){
  ortho_keys[OG_ortho_key$gene[i]]<-OG_ortho_key$orthogroup[i]
}

#print("and made the orthokeys hash")



## read in files
for (sp in species){
  sp_file = paste(aminoacid_dir,"/",orthogroup,"/",sp,"_",orthogroup,"_AA_kmers.csv",sep='')
  unexpressed_name<-paste(sp,"_unexpressed_transcripts",sep='')
  unexpressed_transcripts<-get(unexpressed_name)
  species_OG_AA_kmers <- read.csv(sp_file)
  species_OG_AA_kmers <-species_OG_AA_kmers[!(hash::values(ortho_keys,species_OG_AA_kmers[,1]) %in% unexpressed_transcripts),]
  object_name<-paste(sp,"_OG_AA_Kmers",sep='')
  assign(object_name,species_OG_AA_kmers)
}


#print(c("Test9"))
## calculate change in kmer diversity for each species
for (sp in species){
  object_name<-paste(sp,"_OG_AA_Kmers",sep='')
  sp_OG_AA_kmers <-get(object_name)
  print(ncol(sp_OG_AA_kmers))
  sp_change_kmer_div<-Change_in_kmer_diversity(sp_OG_AA_kmers)
  colnames(sp_change_kmer_div)<-c("seq","Total_H","New_H","H_change" )
  change_kmer_name<-paste(sp,"_change_kmer_div",sep='')
  assign(change_kmer_name,sp_change_kmer_div)
}
#print(ncol(species1_OG_AA_kmers))
#species1_change_kmer_div<-Change_in_kmer_diversity(species1_OG_AA_kmers)
#print("So we pass Species 1")


#print(1)


OG_ortho_classes_file <- paste(orthogroups_info_dir,"/",orthogroup,"/",orthogroup,"_classes.csv",sep='')
OG_ortho_classes<-read.csv(OG_ortho_classes_file)

#print("we read the orthoclasses too")
#####
for (i in 1:nrow(species_pairs)){
  possible_string1 = paste(species_pairs[i,1],species_pairs[i,2],sep=".")
  possible_string2 = paste(species_pairs[i,2],species_pairs[i,1],sep=".")
  
  ## select only the relevant two species colums
  if (possible_string1 %in% colnames(OG_ortho_classes)){
    OG_classes<-as.data.frame(cbind(OG_ortho_classes$paralog,OG_ortho_classes[possible_string1]))
  } else if (possible_string2 %in% colnames(OG_ortho_classes)){
    OG_classes<-as.data.frame(cbind(OG_ortho_classes$paralog,OG_ortho_classes[possible_string2]))
  }
  
  colnames(OG_classes)<-c("OG_orthos","OG_classes")
  OG_classes<-na.omit(OG_classes)
  
  OG_classes_name<-paste(species_pairs[i,1],"_",species_pairs[i,2],"_OG_classes",sep='')
  assign(OG_classes_name,OG_classes)
}
######
#print("made the classes table")

######
for (i in 1:nrow(species_pairs)){
  #print(i)
  OG_classes_name<-paste(species_pairs[i,1],"_",species_pairs[i,2],"_OG_classes",sep='')
  print(OG_classes_name)
  ortho_classes_name<-paste(species_pairs[i,1],"_",species_pairs[i,2],"_ortho_classes",sep='')
  print(ortho_classes_name)
  OG_classes<-get(OG_classes_name)
  
ortho_classes<-hash()
for (i in 1:length(OG_classes$OG_orthos)){
  ortho_classes[OG_classes$OG_orthos[i]]<-OG_classes$OG_classes[i]
}

assign(ortho_classes_name,ortho_classes)
}
######
#print("made the classes hash")



###
for (i in 1:nrow(species_pairs)){
  #print(1)
  species1<-species_pairs[i,1]
  species2<-species_pairs[i,2]
  species1_change_kmer_div_name<-paste(species1,"_change_kmer_div",sep='')
  species1_change_kmer_div<-get(species1_change_kmer_div_name)
  species2_change_kmer_div_name<-paste(species2,"_change_kmer_div",sep='')
  species2_change_kmer_div<-get(species2_change_kmer_div_name)
  ortho_classes_name<-paste(species1,"_",species2,"_ortho_classes",sep='')
  ortho_classes<-get(ortho_classes_name)

#print("1.1")
#print(species1)
#print(species2)
species1_change_kmer_div<-as.data.frame(species1_change_kmer_div)
#print("1.1.1")
species1_change_kmer_div
#print("1.1.1.1")
#print(colnames(species1_change_kmer_div))
#print(is.data.frame(species1_change_kmer_div))
#print(species1_change_kmer_div[["seq"]])
#print("1.1.1.1.1")
species1_gene_classes<-as.vector(hash::values(ortho_classes,as.vector(hash::values(ortho_keys,species1_change_kmer_div[["seq"]]))))

species2_change_kmer_div<-as.data.frame(species2_change_kmer_div)
#print(colnames(species2_change_kmer_div))
species2_gene_classes<-as.vector(hash::values(ortho_classes,as.vector(hash::values(ortho_keys,species2_change_kmer_div$seq))))

#print("made gene class vectors")
  
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

#print(2)
## our seq_list_3 is OG_ortho_key$gene
AA_div_orthogroups<-hash::values(ortho_keys, All_species_change_kmer_div$seq)
#print(3)
All_species_change_kmer_div<-cbind(All_species_change_kmer_div,AA_div_orthogroups)
#print(4)
change_kmer_div_name<-paste(species1,"_",species2,"_all_change_kmer_div",sep="")
#print(5)
assign(change_kmer_div_name,All_species_change_kmer_div)



p<-ggplot(All_species_change_kmer_div, aes(x=AA_div_orthogroups,y=change_perc,fill=species_labels))+
geom_bar(stat="identity", position = position_dodge()) + theme_cowplot() + ylab("Contribution to AA community divergence (%)")+
  theme(axis.text.x = element_text(angle = 90)) + scale_fill_discrete(labels=c(species1, species2))
p

ggfile_name<-paste(output_dir,"/",orthogroup,"/",species1,"_",species2,"_AA_divergence_by_gene.svg",sep = '')
print(ggfile_name)
ggsave(ggfile_name,plot = last_plot(), device="svg")

##plot_file<-paste("new_figures/",plot_name,".tiff",sep='')
##print(plot_file)
##tiff(plot_file)
print(p)
dev.off
}
###


####
for (i in 1:nrow(species_pairs)){
  species1<-species_pairs[i,1]
  species2<-species_pairs[i,2]
  species1_OG_AA_kmers_name<-paste(species1,"_OG_AA_Kmers",sep='')
  species1_OG_AA_kmers<-get(species1_OG_AA_kmers_name)
  species2_OG_AA_kmers_name<-paste(species2,"_OG_AA_Kmers",sep='')
  species2_OG_AA_kmers<-get(species2_OG_AA_kmers_name)
  
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

species1_species2_change_kmer_div_name<-paste(species1,"_",species2,"_all_change_kmer_div",sep='')
species1_species2_all_change_kmer_div<-get(species1_species2_change_kmer_div_name)

OG_cats<-hash()
for (j in 1:nrow(species1_species2_all_change_kmer_div)){
  OG_cats[species1_species2_all_change_kmer_div$seq[j]]<-species1_species2_all_change_kmer_div$gene_classes[j]
}

species1_class_assignments<-hash::values(OG_cats,species1_OG_AA_kmers$seqid)

species2_class_assignments<-hash::values(OG_cats,species2_OG_AA_kmers$seqid)

species1_reorg_kmers<-cbind(species1_ortho_assignments,species1_class_assignments,
                         species1_OG_AA_kmers)
colnames(species1_reorg_kmers)<-c("species_ortho_assignments","species_class_assignments","seqid","kmer")
species1_reorg_kmers_name<-paste(species1,"_reorg_kmers",sep='')
assign(species1_reorg_kmers_name,species1_reorg_kmers)

species2_reorg_kmers<-cbind(species2_ortho_assignments,species2_class_assignments,
                         species2_OG_AA_kmers)
colnames(species2_reorg_kmers)<-c("species_ortho_assignments","species_class_assignments","seqid","kmer")
species2_reorg_kmers_name<-paste(species2,"_reorg_kmers",sep='')
assign(species2_reorg_kmers_name,species2_reorg_kmers)


#print("BUT BAD DOWN HERE?")

print(unique(species1_reorg_kmers$species_class_assignments))
print(unique(species1_reorg_kmers$seqid))
print(unique(species2_reorg_kmers$species_class_assignments))
print(unique(species2_reorg_kmers$seqid))


## Overall similarity
species1_species2_similarity<-Calc_kmer_similarity(species1_reorg_kmers,species2_reorg_kmers)
species1_species2_ortho_effects<-Calc_kmer_disimilarity_effect(species1_reorg_kmers,species2_reorg_kmers)

output_table_file = paste(output_dir,"/",orthogroup,"/",species1,"_",species2,"_dissimilarity_effects_table.csv",sep = '')
write.csv(species1_species2_ortho_effects,output_table_file, row.names = FALSE)



#print("3")
## Minus conserved orthologs
species1_paralog_kmers<-subset(species1_reorg_kmers,species1_reorg_kmers$species_class_assignments != "conserved")
species2_paralog_kmers<-subset(species2_reorg_kmers,species2_reorg_kmers$species_class_assignments != "conserved")
if (nrow(species1_paralog_kmers) == 0 | nrow(species2_paralog_kmers) == 0){paralog_similarity<-NA} else{
  paralog_similarity<-Calc_kmer_similarity(species1_paralog_kmers,species2_paralog_kmers)  
}

#print("4")
## Minus paralogs (all kinds)
species1_orthos_kmers<-subset(species1_reorg_kmers,species1_reorg_kmers$species_class_assignments == "conserved")
print("5.1")
species2_orthos_kmers<-subset(species2_reorg_kmers,species2_reorg_kmers$species_class_assignments == "conserved" )
print("5.2")
if (nrow(species1_orthos_kmers) == 0 | nrow(species2_orthos_kmers) == 0){orthos_similarity<-NA} else{
  orthos_similarity<-Calc_kmer_similarity(species1_orthos_kmers,species2_orthos_kmers)
}
#print("5.3")
## Minus in-paralogs
#Scate_noninparalog_kmers<-subset(Scate_reorg_Kmers,Scate_reorg_Kmers$Scate_class_assignments != "in-paralog")
#Sterg_noninparalog_kmers<-subset(Sterg_reorg_Kmers,Sterg_reorg_Kmers$Sterg_class_assignments != "in-paralog")
#in_paralog_similarity<-Calc_kmer_similarity(Scate_noninparalog_kmers,Sterg_noninparalog_kmers)


## Minus lost-paralogs
#Scate_nonlostparalog_kmers<-subset(Scate_reorg_Kmers,Scate_reorg_Kmers$Scate_class_assignments != "ca-paralog")
#Sterg_nonlostparalog_kmers<-subset(Sterg_reorg_Kmers,Sterg_reorg_Kmers$Sterg_class_assignments != "ca-paralog")
#lost_paralog_similarity<-Calc_kmer_similarity(Scate_nonlostparalog_kmers,Sterg_nonlostparalog_kmers)
#print("5")
dissimilarities<-c(species1_species2_similarity[1],paralog_similarity[1],
                   orthos_similarity[1])
#dissimilarities<-c(species1_species2_similarity[1],paralog_similarity[1],
#                   orthos_similarity[1],in_paralog_similarity[1], 
#                   lost_paralog_similarity[1])
#print("6")
groups<-c("Overall","Paralogs","Orthologs")
#groups<-c("Overall","Paralogs","Orthologs","In-paralogs","Out-paralogs")
expr_class<-rep("kmers_only",length(dissimilarities))

#print("7")

#print("kmer group similarities")

kmer_group_similarities<-cbind.data.frame(groups,dissimilarities,expr_class)
print(kmer_group_similarities)
kmer_group_similarities_name<-paste(species1,'_',species2,'_kmer_group_similarities',sep='')
assign(kmer_group_similarities_name,kmer_group_similarities)


#print("8")
ortho_all_similarity<-subset(kmer_group_similarities, kmer_group_similarities$groups == "Orthologs" | kmer_group_similarities$groups == "Overall")

p<-ggplot(data=ortho_all_similarity, aes(x=groups, y=dissimilarities)) +
  geom_bar(stat="identity", fill=c("cyan3")) + labs(x="Gene Class", y="Dissimilarity")
p

ggfile_name<-paste(output_dir,"/",orthogroup,"/",species1,"_",species2,"_AA_dissimilarity_by_class.svg",sep = '')
print(ggfile_name)
ggsave(ggfile_name,plot = last_plot(), device="svg")

##plot_file<-paste("new_figures/",plot_name,".tiff",sep='')
##print(plot_file)
##tiff(plot_file)
print(p)
dev.off


## First, just add expression counts to the kmer table
## First step, I have to generate averages (and actually provide an alternative strategy for if there is only one sample per species). This is implemented in the several if statements below.
## Sample 1
species1_samples<-subset(sample_assignment,sample_assignment$Species==species1)
##print(species1_samples)
species2_samples<-subset(sample_assignment,sample_assignment$Species==species2)
##print(species2_samples)
if ( nrow(species1_samples) == 0){print("Uh-oh. Looks like there are no expression samples for species1")}
if ( nrow(species1_samples) == 1 ){
  sample1 = str_replace_all(species1_samples[1,1], '-', '.')
  ##print(colnames(OG_expression))
  ##print(sample1)
  sample1_expression<-OG_expression[sample1]
}
if ( nrow(species1_samples) > 1 ){
  samples_tab<-OG_expression[,1]
  for ( j in 1:nrow(species1_samples)){
    ##print(colnames(OG_expression))
    ##print(species1_samples[j,1])
    sample1 = str_replace_all(species1_samples[j,1], '-', '.')
    ##print(sample1)
    sample1_expression<-OG_expression[sample1]
    samples_tab<-cbind(samples_tab,sample1_expression)
  }
  ##print(colnames(samples_tab))
  sample1_expression<-rowMeans(samples_tab[2:ncol(samples_tab)])
  sample1_expression<-as.data.frame(sample1_expression)
}

## And again for sample2
if ( nrow(species2_samples) == 0){print("Uh-oh. Looks like there are no expression samples for species2")}
if ( nrow(species2_samples) == 1 ){
  sample2 = str_replace_all(species2_samples[1,1], '-', '.')
  sample2_expression<-OG_expression[sample2]
}
if ( nrow(species2_samples) > 1 ){
  samples_tab<-OG_expression[,1]
  for ( j in 1:nrow(species2_samples)){
    sample2 = str_replace_all( species2_samples[j,1], '-', '.')
    sample2_expression<-OG_expression[sample2]
    samples_tab<-cbind(samples_tab,sample2_expression)
  }
  ##print(colnames(samples_tab))
  sample2_expression<-rowMeans(samples_tab[2:ncol(samples_tab)])
  sample2_expression<-as.data.frame(sample2_expression)
}

## Second, step create hashes
species1_expression_hash<-hash()
for (j in 1:nrow(OG_expression)){
  species1_expression_hash[OG_expression$Orthogroup[j]]<-sample1_expression[[1]][j]
}

species2_expression_hash<-hash()
for (j in 1:nrow(OG_expression)){
  species2_expression_hash[OG_expression$Orthogroup[j]]<-sample2_expression[[1]][j]
}

#print(6)

species1_reorg_kmers_name<-paste(species1,"_reorg_kmers",sep='')
species1_reorg_kmers<-get(species1_reorg_kmers_name)
species2_reorg_kmers_name<-paste(species2,"_reorg_kmers",sep='')
species2_reorg_kmers<-get(species2_reorg_kmers_name)
## Then generate the expression vectors
species1_expression_count<-hash::values(species1_expression_hash,species1_reorg_kmers$species_ortho_assignments)
species2_expression_count<-hash::values(species2_expression_hash,species2_reorg_kmers$species_ortho_assignments)

#print(6.5)

species1_kmer_expr<-cbind(species1_reorg_kmers,species1_expression_count)
colnames(species1_kmer_expr)<-c("species_ortho_assignments","species_class_assignments","seqid",
                             "kmer","expression")
species2_kmer_expr<-cbind(species2_reorg_kmers,species2_expression_count)
colnames(species2_kmer_expr)<-c("species_ortho_assignments","species_class_assignments","seqid",
                             "kmer","expression")

#print(6.55)

species1_kmer_expr_tab<-Make_kmer_expression_table(species1_kmer_expr)
species1_kmer_expr_tab<-na.omit(species1_kmer_expr_tab)

species2_kmer_expr_tab<-Make_kmer_expression_table(species2_kmer_expr)
species2_kmer_expr_tab<-na.omit(species2_kmer_expr_tab)


#print(6.555)

species1_species2_expr_similarity<-Calc_kmer_expr_similarity(species1_kmer_expr_tab,species2_kmer_expr_tab)[1]
species1_species2_expr_dissimilarity_effects<-Calc_kmer_expr_disimilarity_effects(species1_kmer_expr, species2_kmer_expr, species1_kmer_expr_tab,species2_kmer_expr_tab)

output_table_file = paste(output_dir,"/",orthogroup,"/",species1,"_",species2,"_expression_dissimilarity_effects_table.csv",sep = '')
write.csv(species1_species2_expr_dissimilarity_effects,output_table_file, row.names = FALSE)  
  

#print(7)

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


#print(8)

## Minus conserved orthologs
species1_paralog_kmers_expr<-subset(species1_kmer_expr,species1_kmer_expr$species_class_assignments != "conserved")
species2_paralog_kmers_expr<-subset(species2_kmer_expr,species2_kmer_expr$species_class_assignments != "conserved")
species1_paralog_kmers_expr_tab<-Make_kmer_expression_table(species1_paralog_kmers_expr)
species1_paralog_kmers_expr_tab<-na.omit(species1_paralog_kmers_expr_tab)
species2_paralog_kmers_expr_tab<-Make_kmer_expression_table(species2_paralog_kmers_expr)
species2_paralog_kmers_expr_tab<-na.omit(species2_paralog_kmers_expr_tab)
paralog_similarity_expr<-Calc_kmer_expr_similarity(species1_paralog_kmers_expr_tab, species2_paralog_kmers_expr_tab)[1]

#print(9)
## Minus paralogs (all kinds)
species1_orthos_kmers_expr<-subset(species1_kmer_expr,species1_kmer_expr$species_class_assignments == "conserved")
species2_orthos_kmers_expr<-subset(species2_kmer_expr,species2_kmer_expr$species_class_assignments == "conserved")
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

#print(10)

expr_dissimilarities<-c(species1_species2_expr_similarity[1],orthos_similarity_expr[1],
                        paralog_similarity_expr[1])
groups<-c("Overall","Orthologs","Paralogs")
expr_class<-rep("expressed",rep(length(groups)))

expr_group_similarities<-cbind.data.frame(groups,expr_dissimilarities,expr_class)
colnames(expr_group_similarities)<-c("groups","dissimilarities","expr_class")

kmer_group_similarities_name<-paste(species1,'_',species2,'_kmer_group_similarities',sep='')
kmer_group_similarities<-get(kmer_group_similarities_name)
all_group_similarities<-rbind(kmer_group_similarities,expr_group_similarities)


## Here is the orthos and overall plot
ortho_para_expr_similarity<-subset(all_group_similarities, all_group_similarities$groups == "Orthologs" | all_group_similarities$groups == "Overall")

p<-ggplot(data=ortho_para_expr_similarity, aes(x=groups, y=dissimilarities, fill=expr_class)) +
  geom_bar(stat="identity", position=position_dodge()) + labs(x="Gene Class", y="Dissimilarity")




ggfile_name<-paste(output_dir,"/",orthogroup,"/",species1,"_",species2,"_dissimilarity_plot.svg",sep="")
print(ggfile_name)
ggsave(ggfile_name,plot = last_plot(), device="svg")

ggfile_name<-paste(output_dir,"/",orthogroup,"/",species1,"_",species2,"_dissimilarity_plot.png",sep="")
print(ggfile_name)
ggsave(ggfile_name,plot = last_plot(), device="png")

##plot_file<-paste("new_figures/",plot_name,".tiff",sep='')
##print(plot_file)
##tiff(plot_file)
print(p)
dev.off


classifications<-c()
for (j in 1:nrow(ortho_para_expr_similarity)){
  classification = paste(ortho_para_expr_similarity[j,1],ortho_para_expr_similarity[j,3],sep="_")
  classifications<-c(classifications, classification)
}

output_line = t(ortho_para_expr_similarity[,2])
colnames(output_line) = classifications
print(species1)
print(species2)
print(output_line)

output_table_file = paste(output_dir,"/",orthogroup,"/",species1,"_",species2,"_dissimilarity_table.csv",sep = '')
write.csv(output_line,output_table_file, row.names = FALSE)


}
## /great, That will provide us all the dissimilarity information we could possibly need.













library(DESeq2)
## Now for expression information.
OG_expression_file  = paste(expression_dir,"/",orthogroup,"/",orthogroup,"_tpm10k_expression.csv",sep="")
OG_expression = read.csv(OG_expression_file)

OG_counts_file  = paste(expression_dir,"/",orthogroup,"/",orthogroup,"_count_expression.csv",sep="")
OG_counts <- read.csv(OG_counts_file)

## Make table for plotting expression. Species expression will be based on averages
OG_summary_expr<-data.frame(Orthogroup=character(),Species=character(),logTPM10K=numeric(),sd=numeric())
for (sp in species){
  print(sp)
  species_samples<-subset(sample_assignment,sample_assignment$Species==sp)
  if ( nrow(species_samples) == 0){print("Uh-oh. Looks like there are no expression samples for species")}
  if ( nrow(species_samples) == 1 ){
    sample = str_replace_all(species_samples[1,1], '-', '.')
    sample_expression<-log(OG_expression[sample]+1)
    sample_sd<-rep(0,nrow(OG_expression))
    sample_sd<-as.data.frame(sample_sd)
  }
  if ( nrow(species_samples) > 1 ){
    samples_tab<-OG_expression[,1]
    for ( j in 1:nrow(species_samples)){
      sample = str_replace_all(species_samples[j,1], '-', '.')
      sample_expression<-log(OG_expression[sample]+1)
      samples_tab<-cbind(samples_tab,sample_expression)
    }
    print(colnames(samples_tab))
    sample_expression<-rowMeans(samples_tab[2:ncol(samples_tab)])
    sample_expression<-as.data.frame(sample_expression)
    sample_sd<-apply(samples_tab[2:ncol(samples_tab)],1,sd)
    sample_sd<-as.data.frame(sample_sd)
  }
  sp_col<-as.data.frame(rep(sp,nrow(OG_expression)))
  OG_sp_expr<-as.data.frame(cbind(OG_expression[,1],sp_col,sample_expression,sample_sd))
  colnames(OG_sp_expr)<-c('Orthogroup','Species','logTPM10K','sd')
  OG_summary_expr<-rbind(OG_summary_expr,OG_sp_expr)
}


library("ggplot2")
library("cowplot")

p<-ggplot(OG_summary_expr, aes(x=Orthogroup,y=logTPM10K,fill=Species))+
  geom_bar(stat="identity", position = position_dodge()) + 
  geom_errorbar(aes(ymin=logTPM10K-sd, ymax=logTPM10K+sd), width=0.2,position = position_dodge(0.9)) + 
  theme_cowplot() + theme(axis.text.x = element_text(angle = 90))
p


## replace _ in sample names with .
for (i in 1:nrow(sample_assignment)){
  sample = str_replace_all(sample_assignment[i,1], '-', '.')
  sample_assignment[i,1]<-sample
}

row.names(OG_counts)<-OG_counts$Orthogroup
## Just in case, drop columns in counts that are not in samples.
OG_counts<-subset(OG_counts, select=sample_assignment$Sample)

## DESeq Step
for (i in 1:ncol(OG_counts)){
  OG_counts[,i]<-as.integer(OG_counts[,i])
}


## make counts dataset by combining the VG counts with the sample info
dds1 <- DESeqDataSetFromMatrix(countData = OG_counts, colData = sample_assignment, design = ~ Species)
dds1

## associate orthogroup id info with each row
featureData <- data.frame(gene=rownames(OG_counts))
mcols(dds1) <- DataFrame(mcols(dds1), featureData)
mcols(dds1)

## Filter out rows without at least 10 total reads
##keep <- rowSums(counts(dds1)) >= 10
##dds1 <- dds1[keep,]

## Actually run the differential expression analysis

dds1<-estimateSizeFactors(dds1,type="poscounts")
dds1 <- DESeq(dds1)
res <- results(dds1)
overall_res_name<-paste(output_dir,"/",orthogroup,"/",orthogroup,"_DEseq_results.csv",sep = '')
write.csv(as.data.frame(res), overall_res_name)

for (i in 1:nrow(species_pairs)){
  species1<-species_pairs[i,1]
  species2<-species_pairs[i,2]
  
  contrast_list<-c("Species",species1,species2)
  res_tmp <- results(dds1, contrast=contrast_list)
  res_name<-paste(output_dir,"/",orthogroup,"/",orthogroup,"_",species1,"_",species2,"_DEseq_results.csv",sep = '')
  write.csv(as.data.frame(res_tmp),res_name)
}





