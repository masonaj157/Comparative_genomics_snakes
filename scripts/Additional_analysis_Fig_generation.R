
library(argparser)
library(stringr)
library( treeio ) # devtools::install_github("GuangchuangYu/treeio")
library( ggtree ) # devtools::install_github("GuangchuangYu/ggtree")
library( hutan ) # devtools::install_github("caseywdunn/hutan")


library(parallel)
library(foreach)
library(phytools)
library( tidyverse )
library(magrittr)
library(dplyr)
library(hash)
library(phylolm)



parser <- arg_parser("Rscript to extract dissimiliarity information")
# specify our desired options 
# by default ArgumentParser will add an help option 
parser<-add_argument(parser, "--orthogroup", type="character",
                     help="The focal orthogroup")
parser<-add_argument(parser, "--dated_tree", type="character", default = "0_Inputs/species_tree_dated.nwk",
                     help="File containing the dated species tree")
parser<-add_argument(parser, "--samples", type="character", default = "All",
                     help="Either a comma separated list of samples to use for expression calculation, or the value 'All'. Default is 'All' ")
parser<-add_argument(parser, "--sample_assignment", type="character", default = "0_Inputs/sample_assignment.csv",
                     help="File containing the dated species tree")
parser<-add_argument(parser, "--orthogroups_info_dir", type="character", default = "3.6_OG_classifications",
                     help="Directory where orthogroup info and classifications are. default is 3.6_OG_classifications")
parser<-add_argument(parser, "--read_mapping_dir", type="character", default = "1_read_mapping",
                     help="Directory where counts from readmapping and feature counts can be found. Default = 1_read_mapping")
parser<-add_argument(parser, "--output_dir", type="character", default = "4.1_OG_Routputs",
                     help="Directory where to write output files. default is 4.1_OG_Routputs")
args <- parse_args(parser)


orthogroup = args$orthogroup
dated_tree = args$dated_tree
sample_assignment = args$sample_assignment
samples = args$samples
orthogroups_info_dir = args$orthogroups_info_dir
read_mapping_dir = args$read_mapping_dir
output_dir = args$output_dir


dated_tree_file<-dated_tree
sample_assignment_file<-sample_assignment
################################################################################################

##
gene_expression_tab<-function(dated_tree, samples, sample_assignment_dat, read_mapping_dir){
  if (samples != "All"){
    samples_list<-as.vector(strsplit(samples,',')[[1]])
  } else {samples_list<-as.vector(sample_assignment_dat$Sample)}
  
  final_expression_data<-data.frame(Geneid=character(),Expression=numeric())
  for (species in dated_tree$tip.label){
    print("1")
    species_samples<-subset(sample_assignment_dat,sample_assignment_dat$Species==species)
    print(samples_list)
    sub_samples_list = samples_list[samples_list %in% species_samples$Sample]
    sub_samples_list<-str_replace_all(sub_samples_list,'-','.')
    sub_samples_list<-str_replace_all(sub_samples_list,' ','.')
    
    expression_file<-paste(read_mapping_dir,"/",species,"/",species,"_clrTPM10K_venom_expression.csv",sep='')
    expression_data<-read.csv(expression_file)
    
    for (sample in sub_samples_list){
      for (i in range(8,ncol(expression_data))){
        if (grepl(sample, colnames(expression_data[i]), fixed = TRUE)){
          names(expression_data)[names(expression_data) == colnames(expression_data[i])] <- sample
        }}}
    
    Geneid<-expression_data[,2]
    print(colnames(expression_data))
    print(sub_samples_list)
    expression_data<-expression_data[,sub_samples_list]
    print(head(expression_data))
    
    print(class(expression_data))
    print(is.data.frame(expression_data))
    if (is.data.frame(expression_data)){
      print(head(expression_data))
      expression_data<-as.data.frame(rowMeans(expression_data))
      print(head(expression_data))
    }
    else {expression_data<-as.data.frame(expression_data)}
    species_expression_data<-as.data.frame(cbind(Geneid,expression_data))
    colnames(species_expression_data)<-c("Geneid","Expression")
    final_expression_data<-rbind(final_expression_data,species_expression_data)
  }
  #final_expression_data[,2]<-as.numeric(final_expression_data[,2])
  return(final_expression_data)
}
##


clean_genetree_names = function( dated_tree, gene_tree_obj){
    tmp_gene_tree_obj<-gene_tree_obj
    #print(length(tmp_gene_tree_obj@phylo$tip.label))
    for (tip_num in 1:length(tmp_gene_tree_obj@phylo$tip.label)){
      #print(tip_num)
      new_name = str_split(tmp_gene_tree_obj@phylo$tip.label[tip_num],"_",n=2)[[1]][2]
      #print(new_name)
      tmp_gene_tree_obj@phylo$tip.label[tip_num]<-new_name
    }
  for (tip_num in 1:length(tmp_gene_tree_obj@data$clade_label)){
    if (grepl('_', tmp_gene_tree_obj@data$clade_label[tip_num], fixed = TRUE)){
      new_name = str_split(tmp_gene_tree_obj@data$clade_label[tip_num],"_",n=2)[[1]][2]
      #print(new_name)
      tmp_gene_tree_obj@data$clade_label[tip_num]<-new_name
    }
  }
  return(tmp_gene_tree_obj)
}


##
parse_generax_trees = function( tree_text ){
  tree = treeio::read.nhx( tree_text )
  tree@data<-tree@data[order(tree@data$node),]
  
  # Parse clade labels from phylo object into the data frame
  tree@data$clade_label = c( tree@phylo$tip.label, tree@data$S[(length(tree@phylo$tip.label)+1):length(tree@data$S)] )
  
  # Compara trees sometimes have speciation nodes whose descendants have the 
  # same clade name. This isn't biologically possible, so change such nodes 
  # from speciation events to NA so they don't interfere with tree calibration 
  # and are not used for calculating speciation contrasts
  
  n_nodes = nrow( tree@data )
  n_tips = length( tree@phylo$tip.label )
  internal_nodes = ( n_tips + 1 ):n_nodes
  is_speciation = tree@data$D == "N"
  is_speciation[ is.na( is_speciation ) ] = FALSE
  internal_speciation_nodes = tree@data$node[ ( tree@data$node > n_tips ) & is_speciation ]
  
  # Create a vector of clade names of speciation events, with NA for all other nodes
  speciation_names = rep( NA, n_nodes )
  speciation_names[ internal_speciation_nodes ] = tree@data$clade_label[ internal_speciation_nodes ]
  
  # Loop over the internal speciation nodes
  for( i in internal_speciation_nodes ){
    # get descendent internal nodes
    descendants = hutan::descendants( tree@phylo, i )
    descendants = descendants[ descendants %in% internal_nodes ]
    
    descendant_names = speciation_names[ descendants ]
    
    if( speciation_names[i] %in% descendant_names ){
      tree@data$D[i] = NA
    }
  }
  
  
  # Create a human readable Event column
  tree@data$Event = NA
  tree@data$Event[ tree@data$D == "N" ] = "Speciation"
  tree@data$Event[ tree@data$D == "Y" ] = "Duplication"
  tree@data$Event = factor( tree@data$Event, levels=c( "Speciation", "Duplication" ) )
  
  return( tree )
}
##


add_expression_to_tree = function( tree, expression ){
  tree@data %<>% 
    left_join( expression, by = c( "clade_label" = "Geneid" ) )
  return( tree )
}



#' Adjust branch lengths of a phylogenetic tree to make it ultrametric and 
#' to time calibrate the speciation nodes to fixed values
#' 
#' @param nhx A phylogenetic tree as a treeio::treedata object, with speciation 
#' nodes annotated with a value of "N" in @data$D
#' @param calibration_times A dataframe with two columns: age is the age 
#' of a clade, and clade is the name of the clade
#' @param ... Any additional arguments to pass to ape::chronos()
#' @return A treeio::treedata object if successfully calibrated
#' @export

calibrate_tree = function ( nhx, calibration_times, ... ) {
  
  # Create calibration matrix for speciation nodes
  calibration = 
    nhx@data[ !is.tip.nhx( nhx ), ] %>%
    filter( D == "N" ) %>%
    filter( clade_label %in% calibration_times$clade ) %>%
    left_join( calibration_times, c( "clade_label" = "clade" ) ) %>%
    mutate( age.min = age ) %>%
    mutate( age.max = age ) %>% 
    dplyr::select( node, age.min, age.max ) %>%
    mutate( soft.bounds = NA )
  
  tree = try( 
    ape::chronos( nhx@phylo, calibration=calibration, ... ) 
  )
  
  if( "phylo" %in% class( tree ) ){
    class( tree ) = "phylo"
    nhx@phylo = tree
    return( nhx )
  }
  else{
    return( NA )
  }
}



#' Get a boolean vector corresponding to all tips and internal nodes of 
#' a tree, with value TRUE for tips and FALSE for internal nodes
#' 
#' @param nhx A phylogenetic tree as a treeio::treedata object
#' @return A boolean vector
#' @export
is.tip.nhx = function( nhx ) {
  is.tip = rep( FALSE, nrow( nhx@data ) )
  is.tip[ 1:length( nhx@phylo$tip.label ) ] = TRUE
  is.tip
}



#' Record the distance from a tip for each node in a tree in 
#' the data slot
#' 
#' @param nhx A phylogenetic tree as a treeio::treedata
#' @return A treeio::treedata object, with an additional node_age column
#' in the data slot, or the original object if it was not of class treedata
#' @export
store_node_age = function( nhx ) {
  
  if ( class( nhx ) != "treedata" ) {
    return( nhx )
  }
  
  node_age = hutan::distance_from_tip( nhx@phylo )
  
  # make sure the dataframe is ordered by consecutive nodes
  stopifnot( all( nhx@data$node == 1:length( nhx@data$node ) ) )
  
  nhx@data$node_age = node_age 
  
  return( nhx )
}




mod_dN_dS<-function(read_codeml_object,min_subs){
  substitutions_line<-read_codeml_object@data$subs
  sub_counts<-as.numeric(c(lapply(substitutions_line, function(x) if(str_count(x,"") != 0){length(str_split(x,' / ')[[1]])} else{0})))
  read_codeml_object@data <-read_codeml_object@data %>% add_column(sub_counts)
  #data<-read_codeml_object@data
  mod_dN_vs_dS<-c()
  for (i in 1:nrow(read_codeml_object@data)){
    if (read_codeml_object@data$sub_counts[i] > min_subs){
      mod_dN_vs_dS<-c(mod_dN_vs_dS,read_codeml_object@data$dN_vs_dS[i])
    } else{mod_dN_vs_dS<-c(mod_dN_vs_dS,NA)}
  }
  read_codeml_object@data <-read_codeml_object@data %>% add_column(mod_dN_vs_dS)
  return(read_codeml_object)
}



add_parent_info_to_tree<-function(nhx_object){
  Parent_type<-c()
  ## For tip in tree, find parent node
  for (tip_num in 1:length(nhx_object@phylo$tip.label)){
    tip_path<-nodepath(nhx_object@phylo)[[tip_num]]
    parent<-tip_path[length(tip_path)-1]
    parent_type<-as.character(nhx_object@data[nhx_object@data$node == parent,]$Event)
    Parent_type<-c(Parent_type,parent_type)
  }
  Parent_type<-c(Parent_type,rep(NA,nhx_object@phylo$Nnode))
  nhx_object@data <-nhx_object@data %>% add_column(Parent_type)
  return(nhx_object)
}


################################################################################################
################################################################################################
#########   Written loop to make plots per orthogroup. This may be rewritten to do each series of plots per orthogroup
################################################################################################
################################################################################################

## Establish list of orthogroups
orthogroup_list<-c("OG0000005_mod","OG0000016_mod","OG0000589")

## Start loop
for (orthogroup in orthogroup_list){

## Need to use our gene map and gene classes to create a classes list for each species
ortho_keys_file = paste(orthogroups_info_dir,"/",orthogroup,"/",orthogroup,"_genes.csv",sep='')
OG_ortho_key <- read.csv(ortho_keys_file)

## Build hash to associate genes with orthogroups
ortho_keys<-hash()
for (i in 1:nrow(OG_ortho_key)){
  ortho_keys[OG_ortho_key$gene[i]]<-OG_ortho_key$orthogroup[i]
}
  
## Base Prep is done. Next we will read in files, Etc.
  
################################################################################################
################################################################################################
#########                                Read in Data                                  #########
################################################################################################
################################################################################################      
  
### Read in codeml data
## rst file
rstfile_name<-paste("4.2_OG_codeml_nonambiguous_curated/",orthogroup,"/rst",sep='')
rstfile<-system.file(rstfile_name,"rst",package="treeio")
## mlc file, which for us is just the output file  
mlcfile_name<-paste("4.2_OG_codeml_nonambiguous_curated/",orthogroup,"/",orthogroup,"_output",sep='')
mlcfile<-system.file(mlcfile_name,"mlc",package="treeio")
## Read and create object, which is here called ml.  
ml <- read.codeml(rstfile_name, mlcfile_name)
## Use our modify function to change dN/dS cases where we actually don't have data
ml<-mod_dN_dS(ml,4)

### Now we will read in expression data which we keep associated with the gene tree, 
###just in case we want to use it in that context
##dated_tree_file<-"/Users/mason.501/Desktop/Sistrurus_Bothrops_genomes/0_Inputs/species_tree_dated.nwk"
## dated tree to calibrate genetree
dated_tree<-read.newick(dated_tree_file)
## Need to be say which samples are which species
sample_assignment_dat = read.csv(sample_assignment_file)
## Need to modify some gene names because of how different programs replace certain characters
final_expression_data<-gene_expression_tab(dated_tree, samples, sample_assignment_dat, read_mapping_dir)
for (i in 1:length(final_expression_data$Geneid)){
  if (grepl(':cds', final_expression_data$Geneid[i])){
    replacement<-str_replace(final_expression_data$Geneid[i],':cds','')
    final_expression_data$Geneid[i]<-replacement
  }
}


## read in the gene tree NHX
gene_tree = paste(orthogroups_info_dir,"/",orthogroup,"/",orthogroup,"_expr.nhx",sep='')
gene_tree_obj<-parse_generax_trees(gene_tree)
## Add parent data to the tree. As ind the evolutionary event responsible for the parent node
gene_tree_obj<-add_parent_info_to_tree(gene_tree_obj)
## And modify gene tree names for gene tree time calibration
gene_tree_obj<-clean_genetree_names(dated_tree, gene_tree_obj)

## Add expression data to the tree
gene_tree_obj<-add_expression_to_tree( gene_tree_obj, final_expression_data )

## Now we have gene trees with the codeml and expression data associated with them.
## But to be more useful we sould combine those data
################################################################################################   

##Ugh. This seems excessive, but this is how I wrote the code to combine these data
## Basically pull out codeml data
tip_codeml_data<-as.data.frame(ml@data[1:length(ml@phylo$tip.label),,])
row.names(tip_codeml_data)<-ml@phylo$tip.label

## Basically pull out expression data
tip_expression<-as.data.frame(gene_tree_obj@data[1:length(gene_tree_obj@phylo$tip.label),])
## Add orthogroup id's to the data
OGs<-hash::values(ortho_keys,tip_expression$clade_label)
tip_expression<-cbind(tip_expression,OGs)
## Creating another dataframe... seems a bit unnecessary, but I think this one ends up used hereafter
tip_bar<-data.frame(id=tip_expression$clade_label,expression=tip_expression$Expression,species=tip_expression$S)
tip_bar_OG<-hash::values(ortho_keys,tip_bar$id)
tip_bar<-cbind(tip_bar,tip_bar_OG)

## Need to drop tips from expression tree and drop rows from expression and drop rows from other table
tip_expression<-tip_expression[tip_expression$node %in% tip_codeml_data$node,]

## This just reorders codeml and expression data to match. If they don't match then something is weird.
tip_codeml_data<-tip_codeml_data[order(row.names(tip_codeml_data)),]
tip_expression<-tip_expression[order(row.names(tip_expression)),]
if(all(order(row.names(tip_codeml_data)) == order(row.names(tip_expression)))){
  all_tip_data<-cbind(tip_codeml_data,tip_expression[,1:2],tip_expression[,4:ncol(tip_expression)])
} else {"Well something is wrong here"}

## Rename one column so we don't have two 'S' columns
names(all_tip_data)[15] <- 'Species'
## Inefficient way to make a 'lineage' vector
lineage=c()
for (i in 1:length(all_tip_data$Species)){
  if (all_tip_data$Species[i] == "Balternatus" | all_tip_data$Species[i] == "Bfonsecai" | all_tip_data$Species[i] == "Bcotiara"){
    lineage<-c(lineage,"Bothrops")}
  if (all_tip_data$Species[i] == "Smiliarius" | all_tip_data$Species[i] == "Scatenatus" | all_tip_data$Species[i] == "Stergeminus"){
    lineage<-c(lineage,"Sistrurus")}
}
## Add lineage classifications
all_tip_data<-cbind(all_tip_data,lineage)
## Exclude B. alternatus and S. miliarius as they won't be useful for our pairwise comparisons
all_tip_data<-all_tip_data[all_tip_data$Species != "Balternatus" & all_tip_data$Species != "Smiliarius",]

## Finally for codeml comparisons we will exclude paralogs where we did not have sufficient mutations for dNdS
tips_w_subs<-na.omit(all_tip_data)
row.names(tips_w_subs)<-tips_w_subs$clade_label


################################################################################################
################################################################################################
#########               Section 1: Seq (dN/dS) based comparisons                       #########
################################################################################################
################################################################################################  

## Make ggplot of the codeml tree with dN/dS values as branch colors
## This will basically be used as a visual device for the plot and probably appear in supplement
codeml_tree<-ggtree(ml, aes(color=log(mod_dN_vs_dS)),size=1.25,layout = "fan", open.angle = 90) + #geom_tiplab(size=2) + 
  hexpand(0) +
  scale_color_continuous(name='dN_vs_dS', #limits=c(0, 1.5),
                        oob=scales::squish,
                        low='blue', high='red') +
theme_tree(legend.position=c(0.65, .3),bg='transparent',legend.key.size = unit(0.4, 'cm'),plot.margin = unit(c(0,0,0,0), "cm"))

codeml_tree
#tree_pic_file<-paste("~/Sistrurus_Bothrops_genomes/Figures/Figure3/",orthogroup,"_codeml_dNdS_tree.svg",sep='')
#ggsave(tree_pic_file)  
  
  

###############################################################################################
###############################################################################################
###############################################################################################
## At this stage we have codeml data available and associated with specific tips. Now below is the expression data,
## which we can then compare with dN/dS, etc.

## Makes base split violin plot with dots
dNdS_plot<-ggplot(data=tips_w_subs,aes(x=lineage,y=log(mod_dN_vs_dS),group=Species,fill=Species)) + 
  introdataviz::geom_split_violin(aes(fill=Species), show.legend = FALSE, alpha=0.5,width=1.5) + 
  geom_boxplot(position=position_dodge(0.25),width=0.1,alpha=0.5) +
  geom_point(shape=21, size=2.5, position=position_jitterdodge(
    jitter.width = 0.1, dodge.width = 0.25), show.legend = FALSE) +
  scale_fill_manual(values = c("skyblue1","skyblue3","sienna1","sienna3")) + 
  labs(y="log(dN/dS)",xlab=FALSE) + scale_x_discrete(expand=c(2,0),drop=FALSE) +
  theme_classic() + theme(text = element_text(size=15),axis.title.x = element_blank())

## Need to extract the mapping for jittered points for connecting segments
segment_data_raw<-cbind(layer_data(dNdS_plot,i=3),tips_w_subs$OGs,tips_w_subs$lineage,tips_w_subs$Species)

## Inefficient loop to make the segment dataframe
segment_data<-data.frame(x0=numeric(),x1=numeric(),y0=numeric(),y1=numeric(),lineage=character(),OG=character())
for (i in 1:(nrow(segment_data_raw)-1)){
  for (j in (i+1):nrow(segment_data_raw)){
    #if (i==0 | j ==0){print('Ha') else{
    if ((segment_data_raw$`tips_w_subs$OGs`[i] == segment_data_raw$`tips_w_subs$OGs`[j]) & 
        (segment_data_raw$`tips_w_subs$lineage`[i] == segment_data_raw$`tips_w_subs$lineage`[j])){
      print(i)
      segment_data[nrow(segment_data) + 1,] <- c(as.numeric(as.character(segment_data_raw$x[i])),as.numeric(as.character(segment_data_raw$x[j])),
                                                 segment_data_raw$y[i], segment_data_raw$y[j],
                                                 segment_data_raw$`tips_w_subs$lineage`[i], 
                                                 segment_data_raw$`tips_w_subs$OGs`[i])
    }}
}

## Ughh. I wouldn't have this problem if I were better at this
## Make the x and y columns numeric
segment_data[,1]<-as.numeric( segment_data[,1])
segment_data[,2]<-as.numeric( segment_data[,2])
segment_data[,3]<-as.numeric( segment_data[,3])
segment_data[,4]<-as.numeric( segment_data[,4])

## Alrgiht, add lines connecting orthogroups
dNdS_plot<-dNdS_plot+geom_segment(data=segment_data,aes(x=x0,xend=x1,
                                          y=y0,yend=y1,color=lineage),inherit.aes = FALSE,size=1.25,alpha=0.85) +
  scale_color_manual(values=c("steelblue3","sienna2"))
dNdS_plot + annotation_custom(ggplotGrob(codeml_tree), xmin=-1.7, xmax = 0.8, ymin=-9.75, ymax=5)

dNdS_file<-paste("~/Desktop/Sistrurus_Bothrops_genomes/Figures/Figure3/",orthogroup,"_dNdS_new.svg",sep='')
ggsave(dNdS_file)  


## Need to do T-test and pairwise T-test for orthogroups or Wilcox and pairwise wilcox
shapiro.test(tips_w_subs[tips_w_subs$Species=="Bcotiara",]$mod_dN_vs_dS)
shapiro.test(tips_w_subs[tips_w_subs$Species=="Bfonsecai",]$mod_dN_vs_dS)
shapiro.test(tips_w_subs[tips_w_subs$Species=="Scatenatus",]$mod_dN_vs_dS)
shapiro.test(tips_w_subs[tips_w_subs$Species=="Stergeminus",]$mod_dN_vs_dS)
## So, not enough normality to justify t-tests
Bc_Bf_wilcox<-wilcox.test(tips_w_subs[tips_w_subs$Species=="Bcotiara",]$mod_dN_vs_dS,tips_w_subs[tips_w_subs$Species=="Bfonsecai",]$mod_dN_vs_dS,exact=FALSE)
Sc_St_wilcox<-wilcox.test(tips_w_subs[tips_w_subs$Species=="Scatenatus",]$mod_dN_vs_dS,tips_w_subs[tips_w_subs$Species=="Stergeminus",]$mod_dN_vs_dS,exact=FALSE)
Bc_Bf_wilcox
Sc_St_wilcox

## Subset for each species
Bcoti_tips_w_subs<-tips_w_subs[tips_w_subs$Species=="Bcotiara",]
Bfons_tips_w_subs<-tips_w_subs[tips_w_subs$Species=="Bfonsecai",]
Scate_tips_w_subs<-tips_w_subs[tips_w_subs$Species=="Scatenatus",]
Sterg_tips_w_subs<-tips_w_subs[tips_w_subs$Species=="Stergeminus",]
## Subset to orthologs then make sure orders match
Bcoti_tips_w_subs<-Bcoti_tips_w_subs[Bcoti_tips_w_subs$OGs %in% Bfons_tips_w_subs$OGs,]
Bfons_tips_w_subs<-Bfons_tips_w_subs[Bfons_tips_w_subs$OGs %in% Bcoti_tips_w_subs$OGs,]
Bcoti_tips_w_subs<-Bcoti_tips_w_subs[order(Bcoti_tips_w_subs$OGs),]
Bfons_tips_w_subs<-Bfons_tips_w_subs[order(Bfons_tips_w_subs$OGs),]

Scate_tips_w_subs<-Scate_tips_w_subs[Scate_tips_w_subs$OGs %in% Sterg_tips_w_subs$OGs,]
Sterg_tips_w_subs<-Sterg_tips_w_subs[Sterg_tips_w_subs$OGs %in% Scate_tips_w_subs$OGs,]
Scate_tips_w_subs<-Scate_tips_w_subs[order(Scate_tips_w_subs$OGs),]
Sterg_tips_w_subs<-Sterg_tips_w_subs[order(Sterg_tips_w_subs$OGs),]

wilcox.test(Bcoti_tips_w_subs$mod_dN_vs_dS,Bfons_tips_w_subs$mod_dN_vs_dS, paired=TRUE, exact=FALSE)
wilcox.test(Scate_tips_w_subs$mod_dN_vs_dS,Sterg_tips_w_subs$mod_dN_vs_dS, paired=TRUE, exact=FALSE)


################################################################################################
################################################################################################
#########                         Section 2: Expression                                #########
################################################################################################
################################################################################################  

## We are going to do part of the plots first, then generate ortholog histograms because we will already have
## Orthologs as a table for filtering
expr_plot<-ggplot(data=all_tip_data,aes(x=lineage,y=Expression,group=Species,fill=Species)) + 
  introdataviz::geom_split_violin(aes(fill=Species), show.legend = FALSE, alpha=0.5,width=1) + 
  geom_boxplot(position=position_dodge(0.25),width=0.1,alpha=0.5) +
  geom_point(shape=21, size=2.5, position=position_jitterdodge(
    jitter.width = 0.1, dodge.width = 0.25), show.legend = FALSE) +
  scale_fill_manual(values = c("skyblue1","skyblue3","sienna1","sienna3")) + 
  scale_x_discrete(expand=c(2,0),drop=FALSE) +
  labs(y="Expression (clrTPM10K)",xlab=FALSE) +
  theme_classic() + theme(text = element_text(size=15),axis.title.x = element_blank())

segment_data_raw<-cbind(layer_data(expr_plot,i=3),all_tip_data$OGs,all_tip_data$lineage,all_tip_data$Species)

segment_data<-data.frame(x0=numeric(),x1=numeric(),y0=numeric(),y1=numeric(),lineage=character(),OG=character())
for (i in 1:(nrow(segment_data_raw)-1)){
  for (j in (i+1):nrow(segment_data_raw)){
    #if (i==0 | j ==0){print('Ha') else{
    if ((segment_data_raw$`all_tip_data$OGs`[i] == segment_data_raw$`all_tip_data$OGs`[j]) & 
        (segment_data_raw$`all_tip_data$lineage`[i] == segment_data_raw$`all_tip_data$lineage`[j])){
      print(i)
      segment_data[nrow(segment_data) + 1,] <- c(as.numeric(as.character(segment_data_raw$x[i])),as.numeric(as.character(segment_data_raw$x[j])),
                                                 segment_data_raw$y[i], segment_data_raw$y[j],
                                                 segment_data_raw$`all_tip_data$lineage`[i], 
                                                 segment_data_raw$`all_tip_data$OGs`[i])
    }}
}

segment_data[,1]<-as.numeric( segment_data[,1])
segment_data[,2]<-as.numeric( segment_data[,2])
segment_data[,3]<-as.numeric( segment_data[,3])
segment_data[,4]<-as.numeric( segment_data[,4])

## Now we can easily subset differential expression by orthologs
## Read in DESeq stuff
Scat_Sterg_tip_data<-all_tip_data[all_tip_data$Species=="Scatenatus"|all_tip_data$Species=="Stergeminus",]
Bcoti_Bfons_tip_data<-all_tip_data[all_tip_data$Species=="Bcotiara"|all_tip_data$Species=="Bfonsecai",]

Scate_Sterg_DE_file = paste("4.1_OG_Routputs/",orthogroup,"/",orthogroup,"_Scatenatus_Stergeminus_DEseq_results.csv",sep='')
Scate_Sterg_DE<-read.csv(Scate_Sterg_DE_file)
Scate_Sterg_DE<-Scate_Sterg_DE[Scate_Sterg_DE$X %in% Scat_Sterg_tip_data$OGs,]
Scate_Sterg_DE<-Scate_Sterg_DE[Scate_Sterg_DE$X %in% segment_data[segment_data$lineage=="Sistrurus",]$OG,]


Bcoti_Bfons_DE_file = paste("4.1_OG_Routputs/",orthogroup,"/",orthogroup,"_Bcotiara_Bfonsecai_DEseq_results.csv",sep='')
Bcoti_Bfons_DE<-read.csv(Bcoti_Bfons_DE_file)
Bcoti_Bfons_DE<-Bcoti_Bfons_DE[Bcoti_Bfons_DE$X %in% Bcoti_Bfons_tip_data$OGs,]
Bcoti_Bfons_DE<-Bcoti_Bfons_DE[Bcoti_Bfons_DE$X %in% segment_data[segment_data$lineage=="Bothrops",]$OG,]


Sist_hist<-ggplot(data=Scate_Sterg_DE,aes(x=log2FoldChange)) +
  #geom_histogram(bins=10, color="grey2", fill="grey50") + 
  geom_dotplot(aes(fill=X),stackgroups = TRUE,binpositions="all",show.legend = FALSE,dotsize=2) +
  scale_y_continuous(NULL, breaks = NULL,lim=c(0,0.25)) +
  theme_bw() + theme(legend.text = element_text(size = 6))

Both_hist<-ggplot(data=Bcoti_Bfons_DE,aes(x=log2FoldChange)) +
  #geom_histogram(bins=10, color="grey2", fill="grey50") + 
  geom_dotplot(aes(fill=X),stackgroups = TRUE,binpositions="all",show.legend = FALSE,dotsize=2) + 
  scale_y_continuous(NULL, breaks = NULL,lim=c(0,0.25)) +
  theme_bw() + theme(legend.text = element_text(size = 6))


expr_plot<-expr_plot+geom_segment(data=segment_data,aes(x=x0,xend=x1,
                      y=y0,yend=y1,color=lineage),inherit.aes = FALSE,size=1.25,alpha=0.85) +
  scale_color_manual(values=c("steelblue3","sienna2"))

expr_plot + annotation_custom(ggplotGrob(Both_hist), xmin=-1, xmax = 0.5, ymin=2.5, ymax=6.7) + 
  annotation_custom(ggplotGrob(Sist_hist), xmin=2.5, xmax = 4, ymin=2.5, ymax=6.7)

DE_file<-paste("~/Desktop/Sistrurus_Bothrops_genomes/Figures/Figure4/Mod/",orthogroup,"_DE_distrib_new.svg",sep='')
ggsave(DE_file)  


### Actual statistical comparisons

## Need to do T-test and pairwise T-test for orthogroups or Wilcox and pairwise wilcox
shapiro.test(all_tip_data[all_tip_data$Species=="Bcotiara",]$Expression)
shapiro.test(all_tip_data[all_tip_data$Species=="Bfonsecai",]$Expression)
shapiro.test(all_tip_data[all_tip_data$Species=="Scatenatus",]$Expression)
shapiro.test(all_tip_data[all_tip_data$Species=="Stergeminus",]$Expression)
## So, not enough normality to justify t-tests
Bc_Bf_wilcox<-wilcox.test(all_tip_data[all_tip_data$Species=="Bcotiara",]$Expression,all_tip_data[all_tip_data$Species=="Bfonsecai",]$Expression, exact=FALSE)
Sc_St_wilcox<-wilcox.test(all_tip_data[all_tip_data$Species=="Scatenatus",]$Expression,all_tip_data[all_tip_data$Species=="Stergeminus",]$Expression, exact=FALSE)
Bc_Bf_wilcox
Sc_St_wilcox

## Subset for each species
Bcoti_tips<-all_tip_data[all_tip_data$Species=="Bcotiara",]
Bfons_tips<-all_tip_data[all_tip_data$Species=="Bfonsecai",]
Scate_tips<-all_tip_data[all_tip_data$Species=="Scatenatus",]
Sterg_tips<-all_tip_data[all_tip_data$Species=="Stergeminus",]
## Subset to orthologs then make sure orders match
Bcoti_tips<-Bcoti_tips[Bcoti_tips$OGs %in% Bfons_tips_w_subs$OGs,]
Bfons_tips<-Bfons_tips[Bfons_tips$OGs %in% Bcoti_tips_w_subs$OGs,]
Bcoti_tips<-Bcoti_tips[order(Bcoti_tips$OGs),]
Bfons_tips<-Bfons_tips[order(Bfons_tips$OGs),]

Scate_tips<-Scate_tips[Scate_tips$OGs %in% Sterg_tips$OGs,]
Sterg_tips<-Sterg_tips[Sterg_tips$OGs %in% Scate_tips$OGs,]
Scate_tips<-Scate_tips[order(Scate_tips$OGs),]
Sterg_tips<-Sterg_tips[order(Sterg_tips$OGs),]

wilcox.test(Bcoti_tips$Expression,Bfons_tips$Expression, paired=TRUE)
wilcox.test(Scate_tips$Expression,Sterg_tips$Expression, paired=TRUE)

###############################################################################





ggtree(gene_tree_obj@phylo)+geom_tiplab(size=1.5)+hexpand(.15) +geom_facet(panel = "Trait", data = tip_bar, geom = geom_col,aes(x = expression,fill = tip_bar_OG), orientation = 'y', width = .6)+ 
  theme(legend.text = element_text(size = 2)) + guides(shape = guide_legend(override.aes = list(size = 0.01)))
#tree_exp_file<-paste("~/Sistrurus_Bothrops_genomes/Figures/Figure3/",orthogroup,"_expression_tree.svg",sep='')
#ggsave(tree_exp_file)  

reg_lm<-lm(tips_w_subs$dN_vs_dS~tips_w_subs$Expression)
summary(reg_lm)
kappa_phylolm<-phylolm(data=tips_w_subs,log(dN_vs_dS)~Expression,phy=gene_tree_obj@phylo,model=c("lambda"))
summary(kappa_phylolm)

Scat_Sterg_tip_data<-all_tip_data[all_tip_data$Species=="Scatenatus"|all_tip_data$Species=="Stergeminus",]

Bcoti_Bfons_tip_data<-all_tip_data[all_tip_data$Species=="Bcotiara"|all_tip_data$Species=="Bfonsecai",]

library(patchwork)


###############################################################################
## Code editing Dec. 15 2022 - Adding Dissimilarity information
###############################################################################

## Read in Dissimilarity effects for each species pair

Scate_Sterg_Dissimilarity_Expr<-read.csv(paste("4.1_OG_Routputs/",orthogroup,"/Scatenatus_Stergeminus_expression_dissimilarity_effects_table.csv",sep=''),header=T)

Bcoti_Bfons_Dissimilarity_Expr<-read.csv(paste("4.1_OG_Routputs/",orthogroup,"/Bcotiara_Bfonsecai_expression_dissimilarity_effects_table.csv",sep=''),header=T)


ortho_classifications<-read.csv(paste(orthogroups_info_dir,"/",orthogroup,"/",orthogroup,"_classes.csv",sep=""),header=T)

Scate_Sterg_ortho_hash<-hash()
for(i in 1:nrow(ortho_classifications)){
  Scate_Sterg_ortho_hash[ortho_classifications$paralog[i]] <- ortho_classifications$Stergeminus.Scatenatus[i]
}

Bcoti_Bfons_ortho_hash<-hash()
for(i in 1:nrow(ortho_classifications)){
  Bcoti_Bfons_ortho_hash[ortho_classifications$paralog[i]]<-ortho_classifications$Bfonsecai.Bcotiara[i]
}


## Combine with expression fold changes
Scate_Sterg_DE_file = paste("4.1_OG_Routputs/",orthogroup,"/",orthogroup,"_Scatenatus_Stergeminus_DEseq_results.csv",sep='')
Scate_Sterg_DE<-read.csv(Scate_Sterg_DE_file)
names(Scate_Sterg_DE)[1]<-"orthos"
class<-hash::values(Scate_Sterg_ortho_hash,Scate_Sterg_DE$orthos)
Scate_Sterg_DE<-cbind(Scate_Sterg_DE,class)
Scate_Sterg_DE_Dissimilarity_expression<-merge(Scate_Sterg_DE,Scate_Sterg_Dissimilarity_Expr, by = "orthos")

Bcoti_Bfons_DE_file = paste("4.1_OG_Routputs/",orthogroup,"/",orthogroup,"_Bcotiara_Bfonsecai_DEseq_results.csv",sep='')
Bcoti_Bfons_DE<-read.csv(Bcoti_Bfons_DE_file)
names(Bcoti_Bfons_DE)[1]<-"orthos"
class<-hash::values(Bcoti_Bfons_ortho_hash,Bcoti_Bfons_DE$orthos)
Bcoti_Bfons_DE<-cbind(Bcoti_Bfons_DE,class)
Bcoti_Bfons_DE_Dissimilarity_expression<-merge(Bcoti_Bfons_DE,Bcoti_Bfons_Dissimilarity_Expr, by = "orthos")


## Plot differential expression and logfold change
p1<-ggplot(data=Scate_Sterg_DE_Dissimilarity_expression, aes(x=abs(log2FoldChange),y=disimilarity_effects)) + 
  geom_point(size=5,aes(col=orthos,shape=class),show.legend=T) + #geom_smooth( method = "lm", size = 2) +
  theme_classic()  #+ theme(legend.text = element_text(size = 6))


p2<-ggplot(data=Bcoti_Bfons_DE_Dissimilarity_expression, aes(x=abs(log2FoldChange),y=disimilarity_effects)) + 
  geom_point(size=5,aes(col=orthos,shape=class),show.legend=T) + #geom_smooth( method = "lm", size = 2) +
   theme_classic() # + theme(legend.text = element_text(size = 6))


## Plot differential expression and logfold change
p3<-ggplot(data=Scate_Sterg_DE_Dissimilarity_expression[Scate_Sterg_DE_Dissimilarity_expression$class=="conserved",], aes(x=abs(log2FoldChange),y=disimilarity_effects)) + 
  geom_point(size=5,aes(col=orthos,shape=class),show.legend=T) + #geom_smooth( method = "lm", size = 2) +
   theme_classic()  #+ theme(legend.text = element_text(size = 6))


p4<-ggplot(data=Bcoti_Bfons_DE_Dissimilarity_expression[Bcoti_Bfons_DE_Dissimilarity_expression$class=="conserved",], aes(x=abs(log2FoldChange),y=disimilarity_effects)) + 
  geom_point(size=5,aes(col=orthos,shape=class),show.legend=T) + #geom_smooth( method = "lm", size = 2)
   theme_classic() # + theme(legend.text = element_text(size = 6))






## Plot differential expression and logfold change
p5<-ggplot(data=Scate_Sterg_DE_Dissimilarity_expression[!Scate_Sterg_DE_Dissimilarity_expression$class=="conserved",], aes(x=abs(log2FoldChange),y=disimilarity_effects)) + 
  geom_point(size=5,aes(col=orthos,shape=class),show.legend=T) + #geom_smooth( method = "lm", size = 2) +
   theme_classic()  #+ theme(legend.text = element_text(size = 6))

p6<-ggplot(data=Bcoti_Bfons_DE_Dissimilarity_expression[!Bcoti_Bfons_DE_Dissimilarity_expression$class=="conserved",], aes(x=abs(log2FoldChange),y=disimilarity_effects)) + 
  geom_point(size=5,aes(col=orthos,shape=class),show.legend=T) + #geom_smooth( method = "lm", size = 2) +
   theme_classic() # + theme(legend.text = element_text(size = 6))


(p1 | p2) /
(p3 | p4) /
(p5 | p6)


#DE_file<-paste("~/Sistrurus_Bothrops_genomes/Figures/Dissimilarity/",orthogroup,"_Expression_Dissimilarity.svg",sep='')
#ggsave(DE_file)  



##################################################################################################################################
##################################################################################################################################
##################################################################################################################################
## So, I want to look at both dN/dS in comparison to dissimilarlity and relative expression levels. To do that, we are going to add dissimilarity effects to tables from above.

library(multimode)

Scat_Sterg_tip_data<-all_tip_data[all_tip_data$Species=="Scatenatus"|all_tip_data$Species=="Stergeminus",]

Bcoti_Bfons_tip_data<-all_tip_data[all_tip_data$Species=="Bcotiara"|all_tip_data$Species=="Bfonsecai",]

Scate_Sterg_Dissimilarity<-read.csv(paste("4.1_OG_Routputs/",orthogroup,"/Scatenatus_Stergeminus_dissimilarity_effects_table.csv",sep=''),header=T)
Bcoti_Bfons_Dissimilarity<-read.csv(paste("4.1_OG_Routputs/",orthogroup,"/Bcotiara_Bfonsecai_dissimilarity_effects_table.csv",sep=''),header=T)
Scate_Sterg_Dissimilarity_Expr<-read.csv(paste("4.1_OG_Routputs/",orthogroup,"/Scatenatus_Stergeminus_expression_dissimilarity_effects_table.csv",sep=''),header=T)
Bcoti_Bfons_Dissimilarity_Expr<-read.csv(paste("4.1_OG_Routputs/",orthogroup,"/Bcotiara_Bfonsecai_expression_dissimilarity_effects_table.csv",sep=''),header=T)


Sist_diss<-hash()
for (i in 1:nrow(Scate_Sterg_Dissimilarity)){
  Sist_diss[Scate_Sterg_Dissimilarity$orthos[i]]<-Scate_Sterg_Dissimilarity$disimilarity_effects[i]
}
Both_diss<-hash()
for (i in 1:nrow(Bcoti_Bfons_Dissimilarity)){
  Both_diss[Bcoti_Bfons_Dissimilarity$orthos[i]]<-Bcoti_Bfons_Dissimilarity$disimilarity_effects[i]
}

Sist_expr_diss<-hash()
for (i in 1:nrow(Scate_Sterg_Dissimilarity_Expr)){
  Sist_expr_diss[Scate_Sterg_Dissimilarity_Expr$orthos[i]]<-Scate_Sterg_Dissimilarity_Expr$disimilarity_effects[i]
}
Both_expr_diss<-hash()
for (i in 1:nrow(Bcoti_Bfons_Dissimilarity_Expr)){
  Both_expr_diss[Bcoti_Bfons_Dissimilarity_Expr$orthos[i]]<-Bcoti_Bfons_Dissimilarity_Expr$disimilarity_effects[i]
}

Scat_Sterg_tip_data<-Scat_Sterg_tip_data[Scat_Sterg_tip_data$OGs %in% keys(Sist_diss),]

dissimilarity.effect<-hash::values(Sist_diss,Scat_Sterg_tip_data$OGs)
dissimilarity.expr.effect<-hash::values(Sist_expr_diss,Scat_Sterg_tip_data$OGs)
Scat_Sterg_tip_data<-cbind(Scat_Sterg_tip_data,dissimilarity.effect,dissimilarity.expr.effect)


Bcoti_Bfons_tip_data<-Bcoti_Bfons_tip_data[Bcoti_Bfons_tip_data$OGs %in% keys(Both_diss),]

dissimilarity.effect<-hash::values(Both_diss,Bcoti_Bfons_tip_data$OGs)
dissimilarity.expr.effect<-hash::values(Both_expr_diss,Bcoti_Bfons_tip_data$OGs)
Bcoti_Bfons_tip_data<-cbind(Bcoti_Bfons_tip_data,dissimilarity.effect,dissimilarity.expr.effect)

all_tip_data<-rbind(Scat_Sterg_tip_data,Bcoti_Bfons_tip_data)
tip_data_w_subs<-na.omit(all_tip_data)


p1<-ggplot(data=tip_data_w_subs, aes(x=log(mod_dN_vs_dS),y=log(abs(dissimilarity.effect)),fill=Species,shape=Species)) + 
  geom_point(size=5,show.legend=T) +
  scale_shape_manual(values=c(21,23,23,21))+
  scale_fill_manual(values = c("skyblue1","skyblue3","sienna1","sienna3")) + #geom_smooth( method = "lm", size = 2) +
  theme_classic()  #+ theme(legend.text = element_text(size = 6))

p1

dNdS_diss<-paste("~/Sistrurus_Bothrops_genomes/Figures/Dissimilarity/",orthogroup,"_dNdS_Dissimilarity.svg",sep='')
ggsave(dNdS_diss)  

p3<-ggplot(data=all_tip_data, aes(x=Expression,y=log(abs(dissimilarity.expr.effect)),fill=Species,shape=Species)) + 
  geom_point(size=5,show.legend=T) + #geom_smooth( method = "lm", size = 2) +
  scale_shape_manual(values=c(21,23,23,21))+
  scale_fill_manual(values = c("skyblue1","skyblue3","sienna1","sienna3")) + #geom_smooth( method = "lm", size = 2) +
  theme_classic()  #+ theme(legend.text = element_text(size = 6))

p3

expr_exprdiss<-paste("~/Sistrurus_Bothrops_genomes/Figures/Dissimilarity/",orthogroup,"_Total_Expression_Expr_Dissimilarity.svg",sep='')
ggsave(expr_exprdiss)  

p5<-ggplot(data=all_tip_data, aes(x=Expression,y=dissimilarity.effect,fill=Species,shape=Species)) + 
  geom_point(size=5,show.legend=T) + #geom_smooth( method = "lm", size = 2) +
  scale_shape_manual(values=c(21,23,23,21))+
  scale_fill_manual(values = c("skyblue1","skyblue3","sienna1","sienna3")) + 
  geom_smooth(data=Scat_Sterg_tip_data[Scat_Sterg_tip_data$dissimilarity.effect > -0.04,],
              method = "lm", alpha = .15, aes(x=Expression, y=dissimilarity.effect, group=1), 
              col="sienna2", fill="sienna2")+
  geom_smooth(data=Bcoti_Bfons_tip_data, method = "lm", alpha = .15, 
              aes(x=Expression, y=dissimilarity.effect, group=1), 
              col="skyblue2", fill="skyblue2")+#geom_smooth( method = "lm", size = 2) +
  #geom_smooth( method = "lm", size = 2) +
  theme_classic()  #+ theme(legend.text = element_text(size = 6))

p5

expr_diss<-paste("~/Sistrurus_Bothrops_genomes/Figures/Dissimilarity/",orthogroup,"_Total_Expression_Dissimilarity.svg",sep='')
ggsave(expr_diss)  


lineage=c(rep("Sistrurus",nrow(Scat_Sterg_tip_data)),rep("Bothrops",nrow(Bcoti_Bfons_tip_data)))
all_tip_data<-cbind(all_tip_data,lineage)

p7<-ggplot(data=all_tip_data, aes(x=dissimilarity.effect,fill=lineage)) + 
  geom_density(show.legend=T,alpha=0.45) + 
  scale_fill_manual(values = c("skyblue2","sienna2")) +
  #geom_smooth( method = "lm", size = 2) +
  theme_classic()  #+ theme(legend.text = element_text(size = 6))

p7

diss_distrib<-paste("~/Sistrurus_Bothrops_genomes/Figures/Dissimilarity/",orthogroup,"_Dissimilarity_density.svg",sep='')
ggsave(diss_distrib)  

modetest(Scat_Sterg_tip_data[Scat_Sterg_tip_data$dissimilarity.effect > -0.04,]$dissimilarity.effect,mod0=1)

modetest(Bcoti_Bfons_tip_data[Bcoti_Bfons_tip_data$dissimilarity.effect > -0.04,]$dissimilarity.effect,mod0=1)

summary(lm(data=Scat_Sterg_tip_data[Scat_Sterg_tip_data$dissimilarity.effect > -0.04,],dissimilarity.effect~Expression))

summary(lm(data=Bcoti_Bfons_tip_data[Bcoti_Bfons_tip_data$dissimilarity.effect > -0.04,],dissimilarity.effect~Expression))
}
