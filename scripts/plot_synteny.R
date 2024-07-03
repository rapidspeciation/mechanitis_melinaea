#!/usr/bin/env Rscript
library(argparse)
library(gggenomes)
library(gtools)
library(scales)
library(tibble)
library(tidyr)
library(dplyr)
library(stringr)

# Example script for generating ntSynt synteny ribbon plots using gggenomes
#https://github.com/bcgsc/ntSynt/tree/main/visualization_scripts

#Modified by K Nasvall to work with output from BUSCO full_table.tsv

#This part can be ignored
##############################
# Parse the input arguments
parser <- ArgumentParser(description = "Plot the ntSynt synteny blocks using gggenomes")
parser$add_argument("-s", "--sequences", help = "Input sequence lengths TSV", required = TRUE)
parser$add_argument("-l", "--links", help = "Synteny block links", required = TRUE)
parser$add_argument("--scale", help = "Length of scale bar in bases (default 1 Gbp)", default = 1e9,
                    required = FALSE, type = "double")
parser$add_argument("-p", "--prefix",
                    help = "Output prefix for PNG image (default synteny_gggenomes_plot)", required = FALSE,
                    default = "synteny_gggenomes_plot")

args <- parser$parse_args()
##############################

######################
#prepare input and select taxa
#
#read in renamed busco output table (prepared with with filter_busco_output.R)
links_busco_prel <- read.table(file = paste("../intermediate/", prefix, "_busco_full_table_formated.table", sep = ""))


#switch orientation of whole scaffold if the strand mode is negative rel the busco reference
chain_1 <- links_busco_prel
chain_1$start_init <- chain_1$start
chain_1$end_init <- chain_1$end

chain_1$query <- chain_1$seq_id


Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}


for (i in unique(chain_1$query)) {
  if (Mode(chain_1[chain_1$query==i, "strand"])=="-") {
    print(i)
    
    chain_1[chain_1$query==i, "end"] <- chain_1[chain_1$query==i, "length"] - chain_1[chain_1$query==i, "start"]
    
    chain_1[chain_1$query==i, "start"] <- chain_1[chain_1$query==i, "length"] - chain_1[chain_1$query==i, "end_init"]
    
    chain_1[chain_1$query==i & chain_1$strand=="-", "strand"] <- "temp"
    chain_1[chain_1$query==i & chain_1$strand=="+", "strand"] <- "-"
    chain_1[chain_1$query==i & chain_1$strand=="temp", "strand"] <- "+"
  }
}

links_busco_prel <- chain_1[,c(1:6)]

#remove na
links_busco_prel <- 
  links_busco_prel %>%
  na.omit()

#the order of taxa in the multisynteny plot
unique(links_busco_prel$bin_id)
bin_id_order <- c("ilMelCinx1", "ilDanPlex4", "ilMelLudo1", "ilMelIsoc1", "ilMelMoth8", "ilMelMeno1", "ilMelMars1", "ilMecLysi1", "ilMecMaza1", "ilMecMena1","ilMecMess1", "ilMecPoly1")

##filter away single genes 
#and only include taxa that are in the tree (bin_id_order) 
link_remove <- links_busco_prel[,c("bin_id", "block_id","seq_id")]

#add seq name of ref to each marker in df
marker <- unique(link_remove[link_remove$bin_id==bin_id_order[1], c("seq_id","block_id")])
link_remove <- left_join(link_remove, marker, by="block_id")

#remove markers only present once in a seq (single gene translocations)
link_remove <- 
  link_remove %>%
  group_by(bin_id, seq_id.x, seq_id.y) %>%
  filter(n()>1) %>%
  ungroup()

colnames(link_remove) <- c("bin_id", "block_id", "seq_id", "marker")

links_busco_prel <- inner_join(links_busco_prel, link_remove[,1:3])

#save intermediate
write.table(links_busco_prel, paste("../intermediate/", prefix, "_links_busco_prel_filtered_renamed.table", sep = ""))


#if changing order start here
#links_busco_prel <- read.table("links_busco_prel_filtered_renamed.table")


##order the df
links_temp <- links_busco_prel

#keep only genes found in all taxa
#check number of genes
links_temp[links_temp$block_id %in% Reduce(intersect, split(links_temp$block_id, links_temp$bin_id)),] %>%
  group_by(bin_id) %>%
  count()

#filter
links_temp <- links_temp[links_temp$block_id %in% Reduce(intersect, split(links_temp$block_id, links_temp$bin_id)),]
  
#order seq_id after chromosome number
links_temp <- links_temp[mixedorder(links_temp$seq_id), ]


#order the bin_id after the order in the tree (or list)
links_temp <- links_temp[order(match(links_temp$bin_id,bin_id_order)),]

#order all the markers after the order in the first appearing taxa (ref)
links_temp <- 
links_temp %>%
  mutate(bin_id, bin_id=as.factor(bin_id)) %>%
  arrange(bin_id, forcats::fct_inorder(block_id))

links_temp <- links_temp[order(match(links_temp$bin_id,bin_id_order)),]

#add the columns for the query sequences (all except the first taxa (ref))
sec_taxa <- links_temp[!links_temp$bin_id==bin_id_order[1],]

colnames(sec_taxa) <- paste(colnames(sec_taxa), "2", sep = "")

#add the ref (all except the last taxa), and the query
links_temp <- cbind(links_temp[!links_temp$bin_id==bin_id_order[length(bin_id_order)],], sec_taxa)

#add the sequences that the links should be coloured by
colour_df <-  links_temp[links_temp$bin_id==bin_id_order[1], c("block_id","seq_id")]

colnames(colour_df) <- c("block_id", "colour_block")

links_temp <- left_join(links_temp, colour_df)

#add the strand column, if the strand is opposite in the pairwise comp set strand "-", if they are the same "+"
links_temp$strand_temp <- "+"
links_temp[which(links_temp$strand=="-" & links_temp$strand2=="+"), "strand_temp"] <- "-"
links_temp[which(links_temp$strand=="+" & links_temp$strand2=="-"), "strand_temp"] <- "-"

#use the original strand for column one seq_1 if seq orientation "+", or use links_temp$strand
links_temp$block_ori <- "+"
links_temp$strand <- links_temp$strand_temp

link_busco <- select(links_temp, c("block_id","seq_id","bin_id", "start","end","seq_id2","bin_id2","start2","end2","strand","block_ori","colour_block"))                   

write.table(link_busco, file = paste("../intermediate/", prefix,"_input_full_table_sel_taxa.table", sep = ""))

######################

#if ordering by chr number
# https://stackoverflow.com/questions/32378108/using-gtoolsmixedsort-or-alternatives-with-dplyrarrange
#mixedrank <- function(x) order(gtools::mixedorder(x))
#sequences <- sequences %>%
#  dplyr::arrange(mixedrank(seq_id))

######################
#order by rearrangements, same as the link df
#order reference taxa after chr number and the rest after the order in which they appear in the df
seq_order <- rbind(unique(link_busco[link_busco$bin_id==bin_id_order[1], c("bin_id","seq_id")]), 
                   setNames(unique(link_busco[,c("bin_id2", "seq_id2")]), names(link_busco[,c(3,2)])))

#add length
sequences <- left_join(seq_order, sequences)


write.table(sequences, file = paste("../intermediate/", prefix, "_sequence_order.table", sep = ""))

#for manual fine tuning use the second script ..._fine_tuning.R

######################

#if only plotting and files not in env
#link_busco <- read.table(file = paste(prefix,"input_full_table_sel_taxa.table", sep = "_"))
#sequences <- write.table(sequences, file = paste(prefix, "sequence_order.table", sep = "_"))


###PLOTTING
# Prepare scale bar data frame
scale_bar <- tibble(x = c(0), xend = c(1e9),
                    y = c(0), yend = c(0))

# Infer best units for scale bar
label <- paste(1e9/1e6, "Mb", sep = " ")

label <- paste(args$scale, "bp", sep = " ")
if (args$scale %% 1e9 == 0) {
  label <- paste(args$scale / 1e9, "Gbp", sep = " ")
} else if (args$scale %% 1e6 == 0) {
  label <- paste(args$scale / 1e6, "Mbp", sep = " ")
} else if (args$scale %% 1e3 == 0) {
  label <- paste(args$scale / 1e3, "kbp", sep = " ")
}

#set colours, optional
synt_col <- rev(as.vector(c(pals::kelly(),pals::polychrome()))[2:c(1+length(unique(link_busco$colour_block)))])
#change grey
synt_col[9] <-"#85660D"
synt_col[10] <- "#66B0FF"  
synt_col[11] <- "#f4ed00"  

# Make the ribbon plot - these layers can be fully customized as needed!
make_plot <- function(links, sequences, add_scale_bar = FALSE) {

  p <-  gggenomes(seqs = sequences, links = links)
  plot <- p + theme_gggenomes_clean(base_size = 15) +
    geom_link_line(aes(colour = colour_block), #colour="transparent", offset = 0, 
                   alpha = 0.6, linewidth=0.2) +
    geom_seq(size = 2.8, colour = "grey40") + # draw contig/chromosome lines
    geom_bin_label(size = 4, hjust = 1.5, 
                   #expand_left = 0.5, 
                   fontface = "italic") + # label each bin
    geom_seq_label(aes(label = sub(".*_", "", sequences$seq_id)), vjust = -1.5, hjust = -0.15, size = 2, colour="goldenrod2", fontface = "bold") + # Can add seq labels if desired
    theme(axis.text.x = element_text(size = 25),
          legend.position = "none") +
    scale_fill_manual(values = synt_col,
                      breaks = unique(links$seq_id)
                      ) +
    scale_colour_manual(values = synt_col,
                        breaks = unique(links$seq_id)
    ) +
    guides(fill = guide_legend(ncol = 10),
           colour = guide_legend(title = ""))

  if (add_scale_bar) {
    plot <- plot + geom_segment(data = scale_bar, aes(x = x, xend = xend, y = y, yend = yend),
                                linewidth = 1.5) +
      geom_text(data = scale_bar, aes(x = x + (xend / 2), y = y - 0.3, label = label)) + ggthemr::no_x_axis()
  }

  return(plot)

}


synteny_plot <- make_plot(link_busco, sequences, add_scale_bar = F)
synteny_plot

# Save the ribbon plot
#ggsave(paste(args$prefix, ".png", sep = ""), synteny_plot,
#      units = "cm", widt = 50, height = 20, bg = "white")

ggsave(paste("../output/", prefix, "synteny_gggenomes_prel", ".png", sep = ""), synteny_plot,
       units = "cm", width = 50, height = 50, bg = "white")


