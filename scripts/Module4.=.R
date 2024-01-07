


library(ggplot2)
library(ggridges)
library(hash)
library(stringr)


parser <- arg_parser("Rscript to extract dissimiliarity information")
# specify our desired options 
# by default ArgumentParser will add an help option 
parser<-add_argument(parser, "--input_dissimilarity_directory", type="character",
                     default = "4.1_OG_Routputs",
                     help="The focal orthogroup")
parser<-add_argument(parser, "--keyword_orthogroups", type="character",
                     default = "0_Inputs/keyword_orthogroups.csv",
                     help="The focal orthogroup")
parser<-add_argument(parser, "--genomes_dir", type="character",default = "0_Inputs/genomes",
                     help="Directory containing genomes")
parser<-add_argument(parser, "--output_dir", type="character",
                     help="Directory where to write output files. default is 4.1_OG_Routputs")
args <- parse_args(parser)


input_dissimilarity_directory = args$input_dissimilarity_directory
keyword_orthogroups = args$keyword_orthogroups
genomes_dir = args$genomes_dir
output_dir = args$output_dir


################################################################################################

## This is a very ugly function right now. Perhaps one day I will break it up properly.
ggridges_with_arrow_annotations<-function(dataframe,color_hash){
  
  medians<-c()
  high_percentile<-c()
  for (i in unique(dataframe$group)){
    print(i)
    group_tab<-subset(dataframe,dataframe$group == i)
    print(median(group_tab$Dissimilarity))
    group_median<-median(group_tab$Dissimilarity)
    group_percentile<-quantile(group_tab$Dissimilarity, c(.95)) 
    medians<-c(medians,group_median)
    high_percentile<-c(high_percentile,group_percentile)
  }
  
lines_tab<-data.frame(cbind(unique(dataframe$group),as.numeric(medians),as.numeric(high_percentile)))
lines_tab[,2]<-as.numeric(lines_tab[,2])
lines_tab[,3]<-as.numeric(lines_tab[,3])
colnames(lines_tab)<-c("group","medians","percentiles")
    
    
  p<-ggplot(dataframe, aes(x = Dissimilarity, y = group)) + 
    stat_density_ridges(alpha = 0.7) +
    geom_segment(data = lines_tab, aes(x = medians, xend = medians, y = as.numeric(as.factor(group)),
                                        yend = as.numeric(as.factor(group)) + .9),color = "black")+
    geom_segment(data = lines_tab, aes(x = percentiles, xend = percentiles, y = as.numeric(as.factor(group)),
                                       yend = as.numeric(as.factor(group)) + .9),color = "red")+
    theme_bw()+xlim(0,NA)
  
  plot_groups<-sort(unique(dataframe$group))
  
  group_hash<-hash()
  for (i in 1:length(plot_groups)){
    group_hash[[as.character(plot_groups[i])]]<-i
  }
  
  
  dataframe_toxins<-subset(dataframe,dataframe$toxin_nontoxin == "toxin")
  
  arrow.cor1<-data.frame(x1=double(),y1=double(),x2=double(),y2=double(),group_ver=character())
  text1<-data.frame(x1=double(),y1=double())
  lines1<-data.frame(x1=c(),y1=double(),x2=double(),y2=double())
  for (i in 1:nrow(dataframe_toxins)){
    row<-c(dataframe_toxins$Dissimilarity[i], -1, dataframe_toxins$Dissimilarity[i], 0)
    lines1 <- rbind.data.frame(lines1,row)
    colnames(lines1)=c("x1","y1","x2","y2")
    
    text_row<-c(dataframe_toxins$Dissimilarity[i], 0.6)
    text1<-rbind.data.frame(text1,text_row)
    
    arrow_row<-c(dataframe_toxins$Dissimilarity[i],(hash::values(group_hash,keys=dataframe_toxins$group[i])+0.37),dataframe_toxins$Dissimilarity[i],(hash::values(group_hash,keys=dataframe_toxins$group[i])+0.07),dataframe_toxins$group[i])
    arrow.cor1<-rbind.data.frame(arrow.cor1,arrow_row)
    colnames(arrow.cor1)=c("x1","y1","x2","y2","group_ver")
  }
  
  text1<-cbind(dataframe_toxins$toxin_family,text1)
  colnames(text1)=c("Class","x1","y1")
  
  
  
  
  arrow.cor1<-cbind(arrow.cor1,text1$Class)
  colnames(arrow.cor1)=c("x1","y1","x2","y2","group_ver","Class")
  for (i in 1:4){
    arrow.cor1[,i]<-as.numeric(arrow.cor1[,i])
  }
  
  
  ## Sort table descending to ascending
  arrow.cor1<-arrow.cor1[order(arrow.cor1$x1),]
  
  ## ID arrow_clusters
  ordered_arrows<-order(arrow.cor1$x1,decreasing = TRUE)
  sorted_arrows<-arrow.cor1$x1[ordered_arrows]
  x=0
  arrow_clusters<-c()
  for (i in 1:nrow(arrow.cor1)){
    if (i != 1){
      diff = arrow.cor1$x1[i] - arrow.cor1$x1[i-1]
      #print(diff)
      if ( diff < 0.1){arrow_clusters<-c(arrow_clusters,x)}
      else {
        x=x+1
        arrow_clusters<-c(arrow_clusters,x)
      }}
    else {arrow_clusters<-c(arrow_clusters,x)}}
  
  arrow.cor1<-cbind(arrow.cor1,arrow_clusters)
  
  ## Assign arrow heights
  arrow_cluster = arrow.cor1$arrow_clusters[1]
  if (length(table(arrow.cor1$arrow_clusters == arrow_cluster)) == 2){
    arrow_cluster_max = table(arrow.cor1$arrow_clusters == arrow_cluster)[[2]]*0.05} else{
      arrow_cluster_max = table(arrow.cor1$arrow_clusters == arrow_cluster)[[1]]*0.05  
    }
  y_max = c()
  for (i in 1:nrow(arrow.cor1)){
    if (arrow.cor1$arrow_clusters[i] == arrow_cluster){
      arrow_cluster_max = arrow_cluster_max - .05}
    else {
      arrow_cluster = arrow.cor1$arrow_clusters[i]
      #print(arrow_cluster)
      arrow_cluster_max = ((table(arrow.cor1$arrow_clusters == arrow_cluster)[[2]])-.05)*0.05}
    #print(table(arrow.cor1$clusters == cluster)[[2]])
    y_entry = arrow.cor1$y1[i] + arrow_cluster_max
    y_max = c(y_max, y_entry)
  }
  arrow.cor1<-cbind(arrow.cor1,y_max)
  text_max<-y_max + 0.05
  text1<-text1[order(text1$x1),]
  text1<-cbind(text1,text_max)

  
  p<-p + annotate('segment',x=arrow.cor1$x1,y=arrow.cor1$y_max,xend=arrow.cor1$x2,yend=arrow.cor1$y2,
               lineend = 7, linejoin=c('mitre'), size = 4, col = hash::values(color_hash,keys=text1$Class),
               alpha=0.5,arrow = arrow(length = unit(0.1,"inches"))) +
    #scale_fill_manual(values=toxin_colors) +
    geom_text(data=text1,aes(x=x1, y=text_max,label=Class),
              hjust="left", size = 3)
  
  return(p)
}






make_dissimilarities_table<-function(species_orthos,input_dissimilarity_directory,keyword_orthogroups_tab){
  ##Set table
  Dissimiliarities_table <- data.frame(Overall_kmers_only=numeric(),Ortholog_kmers_only=numeric(),Overall_expressed=numeric(),Orthologs_expressed=numeric())
  row_names_list <-c()
  
  
  ##Build table
  for (i in species_orthos){
    input_file = paste(getwd(),"/",input_dissimilarity_directory,"/",i,"/",
                       species1,'_',species2,"_dissimilarity_table.csv",sep='')
    #print(input_file)
    row = read.csv(input_file)
    row_names_list<-c(row_names_list,i)
    Dissimiliarities_table<-rbind(Dissimiliarities_table,row)
  }
  row.names(Dissimiliarities_table) = row_names_list
  
  
  ## make column of classifications
  toxin_family <-c()
  toxin_nontoxin <-c()
  for (i in 1:nrow(Dissimiliarities_table)){
    if (row.names(Dissimiliarities_table)[i] %in% keyword_orthogroups_tab$Orthogroup){
      val <- as.character(hash::values(keyword_orthogroups_dict, simplify=TRUE, keys=row.names(Dissimiliarities_table)[i])[[1]])
      toxin_family<-c(toxin_family,val)
      toxin_nontoxin<-c(toxin_nontoxin, "toxin")
    } else {
      toxin_family<-c(toxin_family,"nontoxin")
      toxin_nontoxin<-c(toxin_nontoxin, "nontoxin")
    }
  }
  
  Dissimiliarities_table<-cbind(Dissimiliarities_table,toxin_nontoxin,toxin_family)
  
  return(Dissimiliarities_table)
}






################################################################################################

##get list of species from genomes dir
species_files<-list.files(genomes_dir)
species<-c()
for (i in 1:length(species_files)){
  sp = str_split(species_files[i], ".fasta")[[1]][1]
  species<-c(species,sp)
}


species_pairs<-data.frame(species1=character(),species2=character())
for (i in 1:(length(species)-1)){
  print(i)
  for (j in (i+1):length(species)){
    print(j)
    species_pair<-c(species[i],species[j])
    species_pairs<-rbind(species_pairs,species_pair)
  }
}
colnames(species_pairs)<-c("species1","species2")




##read keywords
keyword_orthogroups_tab<-read.csv(keyword_orthogroups)

## build keywords dict
keyword_orthogroups_dict<-hash()
for (i in 1:nrow(keyword_orthogroups_tab)){
  keyword_orthogroups_dict[keyword_orthogroups_tab$Orthogroup[i]]<-keyword_orthogroups_tab$keyword[i]
}

## Read dir
orthogroups <-list.files(path = input_dissimilarity_directory)
for (i in 1:nrow(species_pairs)){
  species1<-species_pairs[i,1]
  species2<-species_pairs[i,2]
  pair_orthos<-c()
  for (ortho in orthogroups){
    ortho_dir <- paste(input_dissimilarity_directory,"/",ortho, sep ='')
    files <- list.files(path = ortho_dir)
    dissimilarity_file = paste(species1,'_',species2,'_dissimilarity_table.csv',sep='')
    if (dissimilarity_file %in% files){pair_orthos<-c(pair_orthos,ortho)}
  }
  pair_orthos_name<-paste(species1,"_",species2,"_orthos",sep='')
  assign(pair_orthos_name,pair_orthos)
  }












for (i in 1:nrow(species_pairs)){
  species1<-species_pairs[i,1]
  species2<-species_pairs[i,2]
  pair_orthos_name<-paste(species1,"_",species2,"_orthos",sep='')
  orthos<-get(pair_orthos_name)

  
  Dissimiliarities_table<-make_dissimilarities_table(orthos,input_dissimilarity_directory,keyword_orthogroups_tab)

  if(species1=="Scatenatus" & species2=="Stergeminus"){
    remove_from_Sistrurus<-c("OG0003361","OG0000344","OG0011740","OG0000381","OG0000076","OG0011741","OG0005687","OG0014164","OG0000609","OG0000285")
    Dissimiliarities_table<-Dissimiliarities_table[!rownames(Dissimiliarities_table) %in% remove_from_Sistrurus,]
  }
  if(species1=="Bcotiara" & species2=="Bfonsecai"){
    remove_from_Bothrops<-c("OG0003361","OG0000344","OG0000381","OG0005811","OG0000076","OG0007625","OG0005687","OG0000609","OG0000285")
    Dissimiliarities_table<-Dissimiliarities_table[!rownames(Dissimiliarities_table) %in% remove_from_Bothrops,]
  }
    
#toxin_colors=c("darkgrey","goldenrod2","yellow","deepskyblue1","red1","green3","darkorange1","green3","deeppink1","darkorange1",
#               "#008080","purple","purple","#fffac8","maroon1","darkgoldenrod","grey","#CC33FF", "palegreen4","mistyrose1","blueviolet","grey2","pink1","khaki")
#toxin_color_names=c("BPP","MYO","SVSP","SVMP","PLA2","CTL","VEGF","CTL","KAZ","VEGF","HYAL","LAAO","LAO","NGF","NUC","PDE","Waprin","CRISP","Vespryn","VF","KUN","Ficolin","PLB","LIPA")
#
toxin_colors<-c("#000000","#CC33FF","green3","#008080","peachpuff","purple","purple","orange","#fffac8",
                        "red3","red1","#800000","royalblue","deepskyblue1","blue1","yellow", "darkorange1",
                        "palegreen4", "mistyrose1","#aa6e28","maroon1","grey","blueviolet","darkgoldenrod","olivedrab1","deeppink1","purple","deepskyblue2")
toxin_color_names<-c("BPP","CRISP","CTL","HYAL","KUN","LAAO","LAO","MYO","NGF","PLA2","Neurotoxic PLA2","PLB","SVMPI",
                       "SVMPII","SVMPIII","SVSP","VEGF","Vespryn","VF","Ficolin","NUC","Waprin","3FTx","PDE","FusedToxin","KAZ","LIPA","SVMP")

color_hash<-hash()
for (i in 1:length(toxin_colors)){
  color_hash[[toxin_color_names[i]]]<-toxin_colors[i]
}



## Table for overall_v_orthologs_expression
orthologs_v_overall_v_expression<-cbind(c(Dissimiliarities_table[,1],Dissimiliarities_table[,2],Dissimiliarities_table[,3]),c(rep("Sequence divergence and geneic variation",length(Dissimiliarities_table[,1])),rep("Sequence divergence",length(Dissimiliarities_table[,2])),rep("Sequence divergence, geneic variation, and expression",length(Dissimiliarities_table[,3]))),c(as.character(Dissimiliarities_table[,5]),as.character(Dissimiliarities_table[,5]),as.character(Dissimiliarities_table[,5])),c(as.character(Dissimiliarities_table[,6]),as.character(Dissimiliarities_table[,6]),as.character(Dissimiliarities_table[,6])))
orthologs_v_overall_v_expression<-as.data.frame(orthologs_v_overall_v_expression)
colnames(orthologs_v_overall_v_expression)<-c("Dissimilarity","group","toxin_nontoxin","toxin_family")
orthologs_v_overall_v_expression$Dissimilarity<-as.numeric(as.character(orthologs_v_overall_v_expression$Dissimilarity))
orthologs_v_overall_v_expression<-na.omit(orthologs_v_overall_v_expression)

orthologs_v_overall_v_expression_tox<-subset(orthologs_v_overall_v_expression,orthologs_v_overall_v_expression$toxin_nontoxin=="toxin")

orthologs_v_overall_v_expression_multi<-subset(orthologs_v_overall_v_expression,
                                             orthologs_v_overall_v_expression$toxin_family=="PLA2"|
                                             orthologs_v_overall_v_expression$toxin_family=="CTL" |
                                             orthologs_v_overall_v_expression$toxin_family=="SVMP"|
                                             orthologs_v_overall_v_expression$toxin_family=="SVSP")





p3<-ggridges_with_arrow_annotations(orthologs_v_overall_v_expression_tox,color_hash)
p3
file_name = paste("5.1_R_Summary_Outputs/",species1,'_',species2,'_dissimiliarity_ridges_tox.svg',sep='')
ggsave(file_name,plot = last_plot(), device="svg")

p4<-ggridges_with_arrow_annotations(orthologs_v_overall_v_expression_multi,color_hash)
p4
file_name = paste("5.1_R_Summary_Outputs/",species1,'_',species2,'_dissimiliarity_ridges_multi.svg',sep='')
ggsave(file_name,plot = last_plot(), device="svg")
 }




############################ Toxin Expression Plots


remove_from_Sistrurus<-c("OG0003361","OG0000344","OG0011740","OG0000381","OG0000076","OG0011741","OG0005687","OG0014164","OG0000609","OG0000285")
species1<-"Scatenatus"
species2<-"Stergeminus"
pair_orthos_name<-paste(species1,"_",species2,"_orthos",sep='')
orthos<-get(pair_orthos_name)

Sc_St_dissimilarities_table<-make_dissimilarities_table(orthos,input_dissimilarity_directory,keyword_orthogroups_tab)

Sc_St_dissimilarities_table<-Sc_St_dissimilarities_table[!rownames(Sc_St_dissimilarities_table) %in% remove_from_Sistrurus,]
Sc_St_toxin_dissimilarities<-subset(Sc_St_dissimilarities_table,Sc_St_dissimilarities_table$toxin_nontoxin == "toxin")

Sc_St_toxin_dissimilarities_multi<-subset(Sc_St_toxin_dissimilarities,
                                            (Sc_St_toxin_dissimilarities$toxin_family == "CTL" | Sc_St_toxin_dissimilarities$toxin_family == "PLA2" | 
                                               Sc_St_toxin_dissimilarities$toxin_family == "SVMP" | Sc_St_toxin_dissimilarities$toxin_family == "SVSP"))

Sc_St_expr_violin_dat<-data.frame(dissimilarity=c(Sc_St_toxin_dissimilarities$Overall_kmers_only,Sc_St_toxin_dissimilarities$Overall_expressed),
                                  group=c(rep("Fgene",nrow(Sc_St_toxin_dissimilarities)),rep("Fexpr",nrow(Sc_St_toxin_dissimilarities))),
                                  Toxin_family=c(Sc_St_toxin_dissimilarities$toxin_family,Sc_St_toxin_dissimilarities$toxin_family))

Sc_St_multi_expr_violin_dat<-data.frame(dissimilarity=c(Sc_St_toxin_dissimilarities_multi$Overall_kmers_only,Sc_St_toxin_dissimilarities_multi$Overall_expressed),
                                  group=c(rep("Fgene",nrow(Sc_St_toxin_dissimilarities_multi)),rep("Fexpr",nrow(Sc_St_toxin_dissimilarities_multi))),
                                  Toxin_family=c(Sc_St_toxin_dissimilarities_multi$toxin_family,Sc_St_toxin_dissimilarities_multi$toxin_family))


Sc_St_expr_wx<-pairwise.wilcox.test(Sc_St_expr_violin_dat$dissimilarity, Sc_St_expr_violin_dat$group, 
                     p.adjust.method = "bonferroni", alternative="less",
                     paired = TRUE)
Sc_St_expr_wx

Sc_St_multi_expr_wx<-pairwise.wilcox.test(Sc_St_multi_expr_violin_dat$dissimilarity, Sc_St_multi_expr_violin_dat$group, 
                                    p.adjust.method = "bonferroni", alternative="less",
                                    paired = TRUE)
Sc_St_multi_expr_wx

## Will need to estimate effect sizes when we can download coin.
##Sc_St_expr_violin_dat %>% wilcox_effsize(dissimilarity~group, paired=TRUE)


p<- ggplot(Sc_St_expr_violin_dat, aes(x=group, y=dissimilarity)) + 
  geom_violin(trim=TRUE,fill="grey") + theme_bw()
p<-p + stat_summary(fun.data="median_hilow", 
                 geom="crossbar", width=0.15 )
p<-p + geom_jitter(shape=16, position=position_jitter(0.2),color = hash::values(color_hash,keys=Sc_St_expr_violin_dat$Toxin_family),size=3)






remove_from_Bothrops<-c("OG0003361","OG0000344","OG0000381","OG0005811","OG0000076","OG0007625","OG0005687","OG0000609","OG0000285")
species1<-"Bcotiara"
species2<-"Bfonsecai"
pair_orthos_name<-paste(species1,"_",species2,"_orthos",sep='')
orthos<-get(pair_orthos_name)

Bc_Bf_dissimilarities_table<-make_dissimilarities_table(orthos,input_dissimilarity_directory,keyword_orthogroups_tab)

Bc_Bf_dissimilarities_table<-Bc_Bf_dissimilarities_table[!rownames(Bc_Bf_dissimilarities_table) %in% remove_from_Bothrops,]
Bc_Bf_toxin_dissimilarities<-subset(Bc_Bf_dissimilarities_table,Bc_Bf_dissimilarities_table$toxin_nontoxin == "toxin")

Bc_Bf_toxin_dissimilarities_multi<-subset(Bc_Bf_toxin_dissimilarities, 
                                            (Bc_Bf_toxin_dissimilarities$toxin_family == "CTL" | Bc_Bf_toxin_dissimilarities$toxin_family == "PLA2" | 
                                               Bc_Bf_toxin_dissimilarities$toxin_family == "SVMP" | Bc_Bf_toxin_dissimilarities$toxin_family == "SVSP"))

Bc_Bf_expr_violin_dat<-data.frame(dissimilarity=c(Bc_Bf_toxin_dissimilarities$Overall_kmers_only,Bc_Bf_toxin_dissimilarities$Overall_expressed),
                                  group=c(rep("Fgene",nrow(Bc_Bf_toxin_dissimilarities)),rep("Fexpr",nrow(Bc_Bf_toxin_dissimilarities))),
                                  Toxin_family=c(Bc_Bf_toxin_dissimilarities$toxin_family,Bc_Bf_toxin_dissimilarities$toxin_family))

Bc_Bf_multi_expr_violin_dat<-data.frame(dissimilarity=c(Bc_Bf_toxin_dissimilarities_multi$Overall_kmers_only,Bc_Bf_toxin_dissimilarities_multi$Overall_expressed),
                                  group=c(rep("Fgene",nrow(Bc_Bf_toxin_dissimilarities_multi)),rep("Fexpr",nrow(Bc_Bf_toxin_dissimilarities_multi))),
                                  Toxin_family=c(Bc_Bf_toxin_dissimilarities_multi$toxin_family,Bc_Bf_toxin_dissimilarities_multi$toxin_family))


Bc_Bf_expr_wx<-pairwise.wilcox.test(Bc_Bf_expr_violin_dat$dissimilarity, Bc_Bf_expr_violin_dat$group, 
                                    p.adjust.method = "bonferroni", alternative="less",
                                    paired = TRUE)
Bc_Bf_expr_wx

Bc_Bf_multi_expr_wx<-pairwise.wilcox.test(Bc_Bf_multi_expr_violin_dat$dissimilarity, Bc_Bf_multi_expr_violin_dat$group, 
                                    p.adjust.method = "bonferroni", alternative="less",
                                    paired = TRUE)
Bc_Bf_multi_expr_wx


p<- ggplot(Bc_Bf_expr_violin_dat, aes(x=group, y=dissimilarity)) + 
  geom_violin(trim=TRUE,fill="grey") + theme_bw()
p<-p + stat_summary(fun.data="median_hilow", 
                    geom="crossbar", width=0.15 )
p<-p + geom_jitter(shape=16, position=position_jitter(0.2),color = hash::values(color_hash,keys=Bc_Bf_expr_violin_dat$Toxin_family),size=3)


############################


remove_from_Sistrurus<-c("OG0003361","OG0000344","OG0011740","OG0000381","OG0000076","OG0011741","OG0005687","OG0014164","OG0000609","OG0000285")
species1<-"Scatenatus"
species2<-"Stergeminus"
pair_orthos_name<-paste(species1,"_",species2,"_orthos",sep='')
orthos<-get(pair_orthos_name)

Sc_St_dissimilarities_table<-make_dissimilarities_table(orthos,input_dissimilarity_directory,keyword_orthogroups_tab)

Sc_St_dissimilarities_table<-Sc_St_dissimilarities_table[!rownames(Sc_St_dissimilarities_table) %in% remove_from_Sistrurus,]
Sc_St_toxin_dissimilarities<-subset(Sc_St_dissimilarities_table,Sc_St_dissimilarities_table$toxin_nontoxin == "toxin")

Sc_St_toxin_dissimilarities_multi<-subset(Sc_St_toxin_dissimilarities, 
                                            (Sc_St_toxin_dissimilarities$toxin_family == "CTL" | Sc_St_toxin_dissimilarities$toxin_family == "PLA2" | 
                                               Sc_St_toxin_dissimilarities$toxin_family == "SVMP" | Sc_St_toxin_dissimilarities$toxin_family == "SVSP"))

Sc_St_mult_par_violin_dat<-data.frame(dissimilarity=c(Sc_St_toxin_dissimilarities_multi$Orthologs_kmers_only,Sc_St_toxin_dissimilarities_multi$Overall_kmers_only),
                                  group=c(rep("sequence",nrow(Sc_St_toxin_dissimilarities_multi)),rep("paralogy",nrow(Sc_St_toxin_dissimilarities_multi))),
                                  Toxin_family=c(Sc_St_toxin_dissimilarities_multi$toxin_family,Sc_St_toxin_dissimilarities_multi$toxin_family))

Sc_St_par_violin_dat<-data.frame(dissimilarity=c(Sc_St_toxin_dissimilarities$Orthologs_kmers_only,Sc_St_toxin_dissimilarities$Overall_kmers_only),
                                      group=c(rep("sequence",nrow(Sc_St_toxin_dissimilarities)),rep("paralogy",nrow(Sc_St_toxin_dissimilarities))),
                                      Toxin_family=c(Sc_St_toxin_dissimilarities$toxin_family,Sc_St_toxin_dissimilarities$toxin_family))

Sc_St_par_wx<-pairwise.wilcox.test(Sc_St_par_violin_dat$dissimilarity, Sc_St_par_violin_dat$group, 
                                   p.adjust.method = "bonferroni", alternative="less",
                                   paired = TRUE)
Sc_St_par_wx

Sc_St_mult_par_wx<-pairwise.wilcox.test(Sc_St_mult_par_violin_dat$dissimilarity, Sc_St_mult_par_violin_dat$group, 
                                    p.adjust.method = "bonferroni", alternative="less",
                                    paired = TRUE)
Sc_St_mult_par_wx

p<- ggplot(Sc_St_par_violin_dat, aes(x=group, y=dissimilarity)) + 
  geom_violin(trim=TRUE,fill="grey") + theme_bw()
p<-p + stat_summary(fun.data="median_hilow", mult=1, 
                 geom="crossbar", width=0.11 )
p<-p + geom_jitter(shape=16, position=position_jitter(0.2),color = hash::values(color_hash,keys=Sc_St_par_violin_dat$Toxin_family),size=3)





Sc_St_par_violin_dat<-data.frame(dissimilarity=c(Sc_St_toxin_dissimilarities$Orthologs_kmers_only,
                                                  Sc_St_toxin_dissimilarities_multi$Orthologs_kmers_only,Sc_St_toxin_dissimilarities_multi$Overall_kmers_only),
                                  group=c(rep("sequence",nrow(Sc_St_toxin_dissimilarities)),
                                          rep("sequence_multi",nrow(Sc_St_toxin_dissimilarities_multi)),
                                          rep("paralogy",nrow(Sc_St_toxin_dissimilarities_multi))),
                                  Toxin_family=c(Sc_St_toxin_dissimilarities$toxin_family,Sc_St_toxin_dissimilarities_multi$toxin_family,
                                                 Sc_St_toxin_dissimilarities_multi$toxin_family))


p<- ggplot(Sc_St_par_violin_dat, aes(x=group, y=dissimilarity)) + 
  geom_violin(trim=TRUE,fill="grey") + theme_bw()
p<-p + stat_summary(fun.data="median_hilow", mult=1, 
                    geom="crossbar", width=0.11 )
p<-p + geom_jitter(shape=16, position=position_jitter(0.2),color = hash::values(color_hash,keys=Sc_St_par_violin_dat$Toxin_family),size=3)







remove_from_Bothrops<-c("OG0003361","OG0000344","OG0000381","OG0005811","OG0000076","OG0007625","OG0005687","OG0000609","OG0000285")
species1<-"Bcotiara"
species2<-"Bfonsecai"
pair_orthos_name<-paste(species1,"_",species2,"_orthos",sep='')
orthos<-get(pair_orthos_name)

Bc_Bf_dissimilarities_table<-make_dissimilarities_table(orthos,input_dissimilarity_directory,keyword_orthogroups_tab)

Bc_Bf_dissimilarities_table<-Bc_Bf_dissimilarities_table[!rownames(Bc_Bf_dissimilarities_table) %in% remove_from_Bothrops,]
Bc_Bf_toxin_dissimilarities<-subset(Bc_Bf_dissimilarities_table,Bc_Bf_dissimilarities_table$toxin_nontoxin == "toxin")

Bc_Bf_toxin_dissimilarities_multi<-subset(Bc_Bf_toxin_dissimilarities, 
                                            (Bc_Bf_toxin_dissimilarities$toxin_family == "CTL" | Bc_Bf_toxin_dissimilarities$toxin_family == "PLA2" | 
                                               Bc_Bf_toxin_dissimilarities$toxin_family == "SVMP" | Bc_Bf_toxin_dissimilarities$toxin_family == "SVSP"))


Bc_Bf_par_violin_dat<-data.frame(dissimilarity=c(Bc_Bf_toxin_dissimilarities$Orthologs_kmers_only,Bc_Bf_toxin_dissimilarities$Overall_kmers_only),
                                  group=c(rep("sequence",nrow(Bc_Bf_toxin_dissimilarities)),rep("paralogy",nrow(Bc_Bf_toxin_dissimilarities))),
                                  Toxin_family=c(Bc_Bf_toxin_dissimilarities$toxin_family,Bc_Bf_toxin_dissimilarities$toxin_family))

Bc_Bf_mult_par_violin_dat<-data.frame(dissimilarity=c(Bc_Bf_toxin_dissimilarities_multi$Orthologs_kmers_only,Bc_Bf_toxin_dissimilarities_multi$Overall_kmers_only),
                                  group=c(rep("sequence",nrow(Bc_Bf_toxin_dissimilarities_multi)),rep("paralogy",nrow(Bc_Bf_toxin_dissimilarities_multi))),
                                  Toxin_family=c(Bc_Bf_toxin_dissimilarities_multi$toxin_family,Bc_Bf_toxin_dissimilarities_multi$toxin_family))


Bc_Bf_par_wx<-pairwise.wilcox.test(Bc_Bf_par_violin_dat$dissimilarity, Bc_Bf_par_violin_dat$group, 
                                    p.adjust.method = "bonferroni", alternative="less",
                                    paired = TRUE)
Bc_Bf_par_wx
wilcoxsign_test(Bc_Bf_par_violin_dat[Bc_Bf_par_violin_dat$group=="sequence",]$dissimilarity ~ Bc_Bf_par_violin_dat[Bc_Bf_par_violin_dat$group=="paralogy",]$dissimilarity, distribution = "exact",
                       alternative = "less")

Bc_Bf_mult_par_wx<-pairwise.wilcox.test(Bc_Bf_mult_par_violin_dat$dissimilarity, Bc_Bf_mult_par_violin_dat$group, 
                                    p.adjust.method = "bonferroni", alternative="less",
                                    paired = TRUE)
Bc_Bf_mult_par_wx
wt<-wilcoxsign_test(Bc_Bf_mult_par_violin_dat[Bc_Bf_mult_par_violin_dat$group=="sequence",]$dissimilarity ~ Bc_Bf_mult_par_violin_dat[Bc_Bf_mult_par_violin_dat$group=="paralogy",]$dissimilarity, distribution = "exact",
                alternative = "less")


p<- ggplot(Bc_Bf_mult_par_violin_dat, aes(x=group, y=dissimilarity)) + 
  geom_violin(trim=TRUE,fill="grey") + theme_bw()
p<-p + stat_summary(fun.data="median_hilow", mult=1, 
                 geom="crossbar", width=0.15 )
p<-p + geom_jitter(shape=16, position=position_jitter(0.2),color = hash::values(color_hash,keys=Bc_Bf_mult_par_violin_dat$Toxin_family),size=3)






Bc_Bf_par_violin_dat<-data.frame(dissimilarity=c(Bc_Bf_toxin_dissimilarities$Orthologs_kmers_only,
                                                  Bc_Bf_toxin_dissimilarities_multi$Orthologs_kmers_only,Bc_Bf_toxin_dissimilarities_multi$Overall_kmers_only),
                                  group=c(rep("sequence",nrow(Bc_Bf_toxin_dissimilarities)),
                                          rep("sequence_multi",nrow(Bc_Bf_toxin_dissimilarities_multi)),
                                          rep("paralogy",nrow(Bc_Bf_toxin_dissimilarities_multi))),
                                  Toxin_family=c(Bc_Bf_toxin_dissimilarities$toxin_family,Bc_Bf_toxin_dissimilarities_multi$toxin_family,
                                                 Bc_Bf_toxin_dissimilarities_multi$toxin_family))


p<- ggplot(Bc_Bf_par_violin_dat, aes(x=group, y=dissimilarity)) + 
  geom_violin(trim=TRUE,fill="grey") + theme_bw()
p<-p + stat_summary(fun.data="median_hilow", mult=1, 
                    geom="crossbar", width=0.11 )
p<-p + geom_jitter(shape=16, position=position_jitter(0.2),color = hash::values(color_hash,keys=Bc_Bf_par_violin_dat$Toxin_family),size=3)




################################################################################################
################################################################################################
################################################################################################

species1<-"Scatenatus"
species2<-"Stergeminus"
pair_orthos_name<-paste(species1,"_",species2,"_orthos",sep='')
orthos<-get(pair_orthos_name)

remove_from_Sistrurus<-c("OG0003361","OG0000344","OG0011740","OG0000381","OG0000076","OG0011741","OG0005687","OG0014164","OG0000609","OG0000285")


Dissimiliarities_table<-make_dissimilarities_table(orthos,input_dissimilarity_directory,keyword_orthogroups_tab)
Dissimiliarities_table<-subset(Dissimiliarities_table,rownames(Dissimiliarities_table) %in% remove_from_Sistrurus == FALSE)

orthologs_v_overall_v_expression<-cbind(c(Dissimiliarities_table[,1],Dissimiliarities_table[,2],Dissimiliarities_table[,3]),c(rep("Sequence divergence and geneic variation",length(Dissimiliarities_table[,1])),rep("Sequence divergence",length(Dissimiliarities_table[,2])),rep("Sequence divergence, geneic variation, and expression",length(Dissimiliarities_table[,3]))),c(as.character(Dissimiliarities_table[,5]),as.character(Dissimiliarities_table[,5]),as.character(Dissimiliarities_table[,5])),c(as.character(Dissimiliarities_table[,6]),as.character(Dissimiliarities_table[,6]),as.character(Dissimiliarities_table[,6])))
orthologs_v_overall_v_expression<-as.data.frame(orthologs_v_overall_v_expression)
colnames(orthologs_v_overall_v_expression)<-c("Dissimilarity","group","toxin_nontoxin","toxin_family")
orthologs_v_overall_v_expression$Dissimilarity<-as.numeric(as.character(orthologs_v_overall_v_expression$Dissimilarity))
orthologs_v_overall_v_expression<-na.omit(orthologs_v_overall_v_expression)

Sc_St_tox_non_orthologs<-subset(orthologs_v_overall_v_expression,orthologs_v_overall_v_expression$group=="Sequence divergence")

Sc_St_tox_non_paralogs<-subset(orthologs_v_overall_v_expression,orthologs_v_overall_v_expression$group=="Sequence divergence and geneic variation")

Sc_St_tox_non_expression<-subset(orthologs_v_overall_v_expression,orthologs_v_overall_v_expression$group=="Sequence divergence, geneic variation, and expression")


p<- ggplot(Sc_St_tox_non_orthologs, aes(x=toxin_nontoxin, y=Dissimilarity)) + 
  geom_violin(trim=TRUE,fill="grey") + theme_bw() + ggtitle("Sequence Divergence")
p<-p + stat_summary(fun.data="median_hilow", mult=1, 
                    geom="crossbar", width=0.11 )
p
wilcox.test(Sc_St_tox_non_orthologs[Sc_St_tox_non_orthologs$toxin_nontoxin=="nontoxin",]$Dissimilarity,
            Sc_St_tox_non_orthologs[Sc_St_tox_non_orthologs$toxin_nontoxin=="toxin",]$Dissimilarity, alternative="less")
#ggsave("~/Desktop/Sistrurus_Bothrops_genomes/Figures/Figure2/Sc_St_Tox_Non_Seq.svg")

p<- ggplot(Sc_St_tox_non_paralogs, aes(x=toxin_nontoxin, y=Dissimilarity)) + 
  geom_violin(trim=TRUE,fill="grey") + theme_bw() + ggtitle("Sequence Divergence Genic Variation")
p<-p + stat_summary(fun.data="median_hilow", mult=1, 
                    geom="crossbar", width=0.11 )
p
wilcox.test(Sc_St_tox_non_paralogs[Sc_St_tox_non_paralogs$toxin_nontoxin=="nontoxin",]$Dissimilarity,
            Sc_St_tox_non_paralogs[Sc_St_tox_non_paralogs$toxin_nontoxin=="toxin",]$Dissimilarity, alternative="less")
#ggsave("~/Desktop/Sistrurus_Bothrops_genomes/Figures/Figure2/Sc_St_Tox_Non_Genes.svg")

p<- ggplot(Sc_St_tox_non_expression, aes(x=toxin_nontoxin, y=Dissimilarity)) + 
  geom_violin(trim=TRUE,fill="grey") + theme_bw() + ggtitle("Expression")
p<-p + stat_summary(fun.data="median_hilow", mult=1, 
                    geom="crossbar", width=0.11 )
p
wilcox.test(Sc_St_tox_non_expression[Sc_St_tox_non_expression$toxin_nontoxin=="nontoxin",]$Dissimilarity,
            Sc_St_tox_non_expression[Sc_St_tox_non_expression$toxin_nontoxin=="toxin",]$Dissimilarity, alternative="less")
#ggsave("~/Desktop/Sistrurus_Bothrops_genomes/Figures/Figure2/Sc_St_Tox_Non_Expression.svg")




orthologs_v_overall_v_expression_tox<-subset(orthologs_v_overall_v_expression,orthologs_v_overall_v_expression$toxin_nontoxin=="toxin")


ggplot(data=orthologs_v_overall_v_expression_tox,aes(x=group, y=Dissimilarity, group = toxin_family, color=toxin_family)) +
  geom_line(size=2) + geom_point(size = 4.5) + theme_classic()

#ggsave("~/Desktop/Sistrurus_Bothrops_genomes/Figures/Figure3/pieces/Sc_St_lines.svg")

#################

species1<-"Bcotiara"
species2<-"Bfonsecai"
pair_orthos_name<-paste(species1,"_",species2,"_orthos",sep='')
orthos<-get(pair_orthos_name)

remove_from_Bothrops<-c("OG0003361","OG0000344","OG0000381","OG0005811","OG0000076","OG0007625","OG0005687","OG0000609","OG0000285")


Dissimiliarities_table<-make_dissimilarities_table(orthos,input_dissimilarity_directory,keyword_orthogroups_tab)
Dissimiliarities_table<-subset(Dissimiliarities_table,rownames(Dissimiliarities_table) %in% remove_from_Bothrops == FALSE)

orthologs_v_overall_v_expression<-cbind(c(Dissimiliarities_table[,1],Dissimiliarities_table[,2],Dissimiliarities_table[,3]),
                                        c(rep("Sequence divergence and geneic variation",length(Dissimiliarities_table[,1])),
                                          rep("Sequence divergence",length(Dissimiliarities_table[,2])),
                                          rep("Sequence divergence, geneic variation, and expression",
                                              length(Dissimiliarities_table[,3]))),
                                        c(as.character(Dissimiliarities_table[,5]),
                                          as.character(Dissimiliarities_table[,5]),
                                          as.character(Dissimiliarities_table[,5])),
                                        c(as.character(Dissimiliarities_table[,6]),
                                          as.character(Dissimiliarities_table[,6]),
                                          as.character(Dissimiliarities_table[,6])))
orthologs_v_overall_v_expression<-as.data.frame(orthologs_v_overall_v_expression)
colnames(orthologs_v_overall_v_expression)<-c("Dissimilarity","group","toxin_nontoxin","toxin_family")
orthologs_v_overall_v_expression$Dissimilarity<-as.numeric(as.character(orthologs_v_overall_v_expression$Dissimilarity))
orthologs_v_overall_v_expression<-na.omit(orthologs_v_overall_v_expression)

Bc_Bf_tox_non_orthologs<-subset(orthologs_v_overall_v_expression,orthologs_v_overall_v_expression$group=="Sequence divergence")

Bc_Bf_tox_non_paralogs<-subset(orthologs_v_overall_v_expression,orthologs_v_overall_v_expression$group=="Sequence divergence and geneic variation")

Bc_Bf_tox_non_expression<-subset(orthologs_v_overall_v_expression,orthologs_v_overall_v_expression$group=="Sequence divergence, geneic variation, and expression")


p<- ggplot(Bc_Bf_tox_non_orthologs, aes(x=toxin_nontoxin, y=Dissimilarity)) + 
  geom_violin(trim=TRUE,fill="grey") + theme_bw() + ggtitle("Sequence Divergence")
p<-p + stat_summary(fun.data="median_hilow", mult=1, 
                    geom="crossbar", width=0.11 )
p
wilcox.test(Bc_Bf_tox_non_orthologs[Bc_Bf_tox_non_orthologs$toxin_nontoxin=="nontoxin",]$Dissimilarity,
            Bc_Bf_tox_non_orthologs[Bc_Bf_tox_non_orthologs$toxin_nontoxin=="toxin",]$Dissimilarity, alternative="less")
#ggsave("~/Desktop/Sistrurus_Bothrops_genomes/Figures/Figure2/Bc_Bf_Tox_Non_Seq.svg")

p<- ggplot(Bc_Bf_tox_non_paralogs, aes(x=toxin_nontoxin, y=Dissimilarity)) + 
  geom_violin(trim=TRUE,fill="grey") + theme_bw() + ggtitle("Sequence Divergence Genic Variation")
p<-p + stat_summary(fun.data="median_hilow", mult=1, 
                    geom="crossbar", width=0.11 )
p
wilcox.test(Bc_Bf_tox_non_paralogs[Bc_Bf_tox_non_paralogs$toxin_nontoxin=="nontoxin",]$Dissimilarity,
            Bc_Bf_tox_non_paralogs[Bc_Bf_tox_non_paralogs$toxin_nontoxin=="toxin",]$Dissimilarity, alternative="less")
ggsave("~/Desktop/Sistrurus_Bothrops_genomes/Figures/Figure2/Bc_Bf_Tox_Non_Genes.svg")

p<- ggplot(Bc_Bf_tox_non_expression, aes(x=toxin_nontoxin, y=Dissimilarity)) + 
  geom_violin(trim=TRUE,fill="grey") + theme_bw() + ggtitle("Expression")
p<-p + stat_summary(fun.data="median_hilow", mult=1, 
                    geom="crossbar", width=0.11 )
p
wilcox.test(Bc_Bf_tox_non_expression[Bc_Bf_tox_non_expression$toxin_nontoxin=="nontoxin",]$Dissimilarity,
            Bc_Bf_tox_non_expression[Bc_Bf_tox_non_expression$toxin_nontoxin=="toxin",]$Dissimilarity, alternative="less")
ggsave("~/Desktop/Sistrurus_Bothrops_genomes/Figures/Figure2/Bc_Bf_Tox_Non_Expression.svg")



orthologs_v_overall_v_expression_tox<-subset(orthologs_v_overall_v_expression,orthologs_v_overall_v_expression$toxin_nontoxin=="toxin")


ggplot(data=orthologs_v_overall_v_expression_tox,aes(x=group, y=Dissimilarity, group = toxin_family, color=toxin_family)) +
  geom_line(size=2) + geom_point(size = 4.5) + theme_classic()


ggsave("~/Desktop/Sistrurus_Bothrops_genomes/Figures/Figure3/pieces/Bc_Bf_lines.svg")


################################################################################################
################################################################################################
################################################################################################



remove_from_Sistrurus<-c("OG0003361","OG0000344","OG0011740","OG0000381","OG0000076","OG0011741","OG0005687","OG0014164","OG0000609","OG0000285")
species1<-"Scatenatus"
species2<-"Stergeminus"
pair_orthos_name<-paste(species1,"_",species2,"_orthos",sep='')
orthos<-get(pair_orthos_name)

Sc_St_dissimilarities_table<-make_dissimilarities_table(orthos,input_dissimilarity_directory,keyword_orthogroups_tab)

Sc_St_dissimilarities_table<-Sc_St_dissimilarities_table[!rownames(Sc_St_dissimilarities_table) %in% remove_from_Sistrurus,]
Sc_St_toxin_dissimilarities<-subset(Sc_St_dissimilarities_table,Sc_St_dissimilarities_table$toxin_nontoxin == "toxin")

Sc_St_toxin_dissimilarities_multi<-subset(Sc_St_toxin_dissimilarities,Sc_St_toxin_dissimilarities$Overall_kmers_only!=Sc_St_toxin_dissimilarities$Orthologs_kmers_only & 
                                            (Sc_St_toxin_dissimilarities$toxin_family == "CTL" | Sc_St_toxin_dissimilarities$toxin_family == "PLA2" | 
                                               Sc_St_toxin_dissimilarities$toxin_family == "SVMP" | Sc_St_toxin_dissimilarities$toxin_family == "SVSP"))

Sc_St_expr_violin_dat<-data.frame(dissimilarity=c(Sc_St_toxin_dissimilarities_multi$Orthologs_expressed,Sc_St_toxin_dissimilarities_multi$Overall_expressed),
                                  group=c(rep("orthologs",nrow(Sc_St_toxin_dissimilarities_multi)),rep("paralogs",nrow(Sc_St_toxin_dissimilarities_multi))),
                                  Toxin_family=c(Sc_St_toxin_dissimilarities_multi$toxin_family,Sc_St_toxin_dissimilarities_multi$toxin_family))


#Sc_St_expr_wx<-pairwise.wilcox.test(Sc_St_expr_violin_dat$dissimilarity, Sc_St_expr_violin_dat$group, 
#                                    p.adjust.method = "bonferroni", alternative="less",
#                                    paired = TRUE)

p<- ggplot(Sc_St_expr_violin_dat, aes(x=group, y=dissimilarity)) + 
  geom_violin(trim=TRUE,fill="grey") + theme_bw() + ggtitle("Sc St Expression")
p<-p + stat_summary(fun.data="median_hilow", mult=1, 
                    geom="crossbar", width=0.11 )
p<-p + geom_jitter(shape=16, position=position_jitter(0.2),color = hash::values(color_hash,keys=Sc_St_expr_violin_dat$Toxin_family),size=3)
p
ggsave("~/Desktop/Sistrurus_Bothrops_genomes/Figures/Figure2/Sc_St_Single_Multi_Expression.svg")








remove_from_Bothrops<-c("OG0003361","OG0000344","OG0000381","OG0005811","OG0000076","OG0007625","OG0005687","OG0000609","OG0000285")
species1<-"Bcotiara"
species2<-"Bfonsecai"
pair_orthos_name<-paste(species1,"_",species2,"_orthos",sep='')
orthos<-get(pair_orthos_name)

Bc_Bf_dissimilarities_table<-make_dissimilarities_table(orthos,input_dissimilarity_directory,keyword_orthogroups_tab)

Bc_Bf_dissimilarities_table<-Bc_Bf_dissimilarities_table[!rownames(Bc_Bf_dissimilarities_table) %in% remove_from_Bothrops,]
Bc_Bf_toxin_dissimilarities<-subset(Bc_Bf_dissimilarities_table,Bc_Bf_dissimilarities_table$toxin_nontoxin == "toxin")

Bc_Bf_toxin_dissimilarities_multi<-subset(Bc_Bf_toxin_dissimilarities,Bc_Bf_toxin_dissimilarities$Overall_kmers_only!=Bc_Bf_toxin_dissimilarities$Orthologs_kmers_only & 
                                            (Bc_Bf_toxin_dissimilarities$toxin_family == "CTL" | Bc_Bf_toxin_dissimilarities$toxin_family == "PLA2" | 
                                               Bc_Bf_toxin_dissimilarities$toxin_family == "SVMP" | Bc_Bf_toxin_dissimilarities$toxin_family == "SVSP"))

Bc_Bf_expr_violin_dat<-data.frame(dissimilarity=c(Bc_Bf_toxin_dissimilarities_multi$Orthologs_expressed,Bc_Bf_toxin_dissimilarities_multi$Overall_expressed),
                                  group=c(rep("orthologs",nrow(Bc_Bf_toxin_dissimilarities_multi)),rep("paralogs",nrow(Bc_Bf_toxin_dissimilarities_multi))),
                                  Toxin_family=c(Bc_Bf_toxin_dissimilarities_multi$toxin_family,Bc_Bf_toxin_dissimilarities_multi$toxin_family))


Bc_Bf_expr_wx<-pairwise.wilcox.test(Bc_Bf_expr_violin_dat$dissimilarity, Bc_Bf_expr_violin_dat$group, 
                                    p.adjust.method = "bonferroni", alternative="less",
                                    paired = TRUE)

p<- ggplot(Bc_Bf_expr_violin_dat, aes(x=group, y=dissimilarity)) + 
  geom_violin(trim=TRUE,fill="grey") + theme_bw() + ggtitle("Bc Bf Expression")
p<-p + stat_summary(fun.data="median_hilow", mult=1, 
                    geom="crossbar", width=0.15 )
p<-p + geom_jitter(shape=16, position=position_jitter(0.2),color = hash::values(color_hash,keys=Bc_Bf_expr_violin_dat$Toxin_family),size=3)
p
ggsave("~/Desktop/Sistrurus_Bothrops_genomes/Figures/Figure2/Bc_Bf_Single_Multi_Expression.svg")



#########################################################################
#########################################################################
#########################################################################








