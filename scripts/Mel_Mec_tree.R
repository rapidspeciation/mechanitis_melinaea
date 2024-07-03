#script to plot tree used in synteny plot

library(ggplot2)
library(ggpubr)
library(phytools)
library(dplyr)

#if (!require("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install(version = "3.17")
#
#BiocManager::install("ggtree")
library(ggtree)

#read in the tree
tree <- read.newick("../../synteny/input/Mel_Mec_top_ind.tree")

#check tree
ggtree(ladderize(tree), branch.length = "none", right = T) +
  geom_nodelab() +
  geom_tiplab() +
  xlim(c(0,10))

#read in the order of tips in the tree
bin_id_order <- unique(read.table("../intermediate/Mec_Mel_sequence_order_mod6.table")[,1])

# taxa_names_id <- data.frame(cbind(species_name=c("Melitaea cinxia", "Danaus plexippus", "Melinaea ludovica", "Melinaea isocomma", "Melinaea mothone", "Melinaea menophilus", "Melinaea marsaeus", "Mechanitis macrinus", "Mechanitis mazaeus", "Mechanitis menapis", "Mechanitis messenoides", "Mechanitis polymnia"), taxa=bin_id_order))

#make a df with tiplabels and the names to display in the tree. 
#OBS check that the order is correct in the df!!!
taxa_names_id <- data.frame(cbind(species_name=c("Melitaea cinxia ", "Danaus plexippus ", "Mel. ludovica ", "Mel. isocomma ", "Mel. mothone ", "Mel. menophilus ", "Mel. marsaeus ", "Mec. macrinus ", "Mec. mazaeus ", "Mec. menapis ", "Mec. messenoides ", "Mec. polymnia "), taxa=bin_id_order))

new_tiplabel <- taxa_names_id[,c("taxa", "species_name")]
colnames(new_tiplabel) <- c("tiplabel", "taxa_name")

#read in table with chromosome number
nr_chr_per_taxa.df <- read.table(file = "../../../../ithomiini_proj/00_documentation/assembly_stats_mec_mel_plot.csv", header=F, fill = NA, sep = ",")

colnames(nr_chr_per_taxa.df) <- c("tiplabel", "autosomes", "sexchr", "genome_size", "in_chr", "percent", "nr_sex_chr", "nr_z")

nr_chr_per_taxa.df$tot_chr <- nr_chr_per_taxa.df$autosomes + nr_chr_per_taxa.df$nr_z

nr_chr_per_taxa.df

#order the tips in the tree
tree_rotated <- ape::rotateConstr(tree, c(rev(new_tiplabel$tiplabel)))

#plot tree
p_clado <-   
  ggtree(tree_rotated, ladderize = F, branch.length="none", layout="rectangular", right = T) %<+% nr_chr_per_taxa.df %<+% new_tiplabel +
  geom_tiplab(aes(label=paste(paste0('italic("', taxa_name, '")'), "(",tot_chr,")", sep = " ")), size=8, parse=T, offset = 0.4, color = "black")  +
  #geom_nodelab(nudge_x = -0.4,nudge_y = -0.3, size=5) +
  #geom_label(aes(label=Chr_nr), size=5, nudge_x = -0.1, nudge_y = 0,  fill="yellow") +
  xlim(0, 17) +
  theme(plot.margin = unit(c(5,20,5,5), "pt"),
        plot.background = element_rect(fill = "transparent",
                                       colour = NA),
        panel.background = element_rect(fill = "transparent",
                                       colour = NA),
  )


#save tree
ggsave(p_clado,
       bg = "transparent",
       filename = paste("../output/synteny_tree_v3.png", sep = ""),
       device = "png",
       height = 20, width = 7
)


