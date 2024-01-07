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

parser <- arg_parser("Rscript to extract dissimiliarity information")
# specify our desired options 
# by default ArgumentParser will add an help option 
parser<-add_argument(parser, "--orthogroup", type="character",
                     help="The focal orthogroup")
parser<-add_argument(parser, "--species", type="character",
                     help="Comma separated pair of species. Example: Species1,Species2")
parser<-add_argument(parser, "--samples", type="character", 
                     help="Comma separated pair of samples (one from each species) that can be found in the featureCounts output. Example: Sample1,Sample2. Can also use 'Average' for one or both species.")
parser<-add_argument(parser, "--aminoacid_dir", type="character", default = "3.7_OG_AA_kmers",
                     help="Directory where there are amino acid kmer counts for the orthogroup. default is 3.7_OG_AA_kmers")
parser<-add_argument(parser, "--orthogroups_info_dir", type="character", default = "3.6_OG_classifications",
                     help="Directory where orthogroup info and classifications are. default is 3.6_OG_classifications")
parser<-add_argument(parser, "--expression_dir", type="character", default = "3.8_OG_expression",
                     help="Directory where orthogroup info and classifications are. default is 3.8_OG_expression")
parser<-add_argument(parser, "--geneorder_dir", type="character", default = "3.9_OG_geneorder",
                     help="Directory containing a csv file of gene positions. default is 3.9_OG_geneorder")
parser<-add_argument(parser, "--decostar_dir", type="character", default = "3.10_OG_decostar",
                     help="Directory containing a csv file of gene positions. default is 3.10_OG_decostar")
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
decostar_dir = args$decostar_dir
output_dir = args$output_dir


################################################################################################

decostar_gene_order_file = paste(decostar_dir,"/",orthogroup,"/output/ancestral_positions_table.csv", sep='')
decostar_gene_order<-read.csv(decostar_gene_order_file)

#gene_order<-gene_order[-c(14),]

#dummies<-make_alignment_dummies(gene_order, aes(xmin = start, xmax = end, y = molecule, id = gene), on = gene_order$gene[1])

#gene_order<-define_seq_groups(gene_order)

decostar_gene_order$species<-factor(decostar_gene_order$species, levels = c("Smiliarius-Scatenatus-Stergeminus","Scatenatus-Stergeminus","Smiliarius", "Scatenatus", "Stergeminus"))
genes_plot<-ggplot(decostar_gene_order, aes(xmin = start, xmax = end, y = species, fill = gene, forward=orientation)) +
  geom_gene_arrow() +
  #geom_blank(data = dummies) +
  #facet_grid( scales = "free", space = "free")+
  facet_wrap(~ species, scales = "free", ncol = 1) +
  #scale_fill_brewer(palette = "Set3") +
  theme_genes()


plot_file = paste(output_dir,"/",orthogroup,"/",orthogroup,"_genes_plot.png",sep="")
print(genes_plot)
ggsave(plot_file,p, device="png")

################################################################################################





