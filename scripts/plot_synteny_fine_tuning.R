#!/usr/bin/env Rscript
library(gggenomes)
library(gtools)
library(scales)
library(tibble)
library(tidyr)
library(dplyr)
library(gtools)

#Modified by K Nasvall to work with output from BUSCO full_table.tsv
#for manual fine-tuning after creating "sequence" and "link_busco" files with 
#the script plot_synteny_blocks_gggenomes_busco_v1


#variables
prefix <- "Mec_Mel"
#manually change the sequence_order.table and save with new name
#read in the new table
sequences <- read.table(file = paste("../intermediate/", prefix,"_sequence_order_mod6.table", sep = ""))


#get a vector of sequences in order
bin_id_order <- unique(sequences$bin_id)

#read in the link table
link_busco <- read.table(file = paste("../intermediate/",prefix,"_input_full_table_sel_taxa.table", sep = ""))

link_busco$seq_id <- paste(link_busco$bin_id, link_busco$seq_id, sep = "_")
link_busco$seq_id2 <- paste(link_busco$bin_id2, link_busco$seq_id2, sep = "_")
sequences$seq_id <- paste(sequences$bin_id, sequences$seq_id, sep = "_")

########################
#reorient the scaffolds
chain_1 <- link_busco
chain_1$start_init <- chain_1$start
chain_1$end_init <- chain_1$end

chain_1$start_init2 <- chain_1$start2
chain_1$end_init2 <- chain_1$end2

chain_1 <- left_join(chain_1, sequences)
seq2 <- sequences
colnames(seq2) <- c("bin_id2", "seq_id2", "length2")
chain_1 <- left_join(chain_1, seq2)


chain_1$query <- chain_1$seq_id
chain_1$query2 <- chain_1$seq_id2


#reorient the same chr (example the Z) in multiple taxa (x, etc), the number is the appearance in the tree
#paste(bin_id_order[c(5, 10)], "Z", sep = "_") is the name of the chromosome, 
for (i in c(paste(bin_id_order[c(5, 10)], "Z", sep = "_"))) {
  print(i)
  
  chain_1[chain_1$query==i, "end"] <- chain_1[chain_1$query==i, "length"] - chain_1[chain_1$query==i, "start"]
  chain_1[chain_1$query2==i, "end2"] <- chain_1[chain_1$query2==i, "length2"] - chain_1[chain_1$query2==i, "start2"]
  
  chain_1[chain_1$query==i, "start"] <- chain_1[chain_1$query==i, "length"] - chain_1[chain_1$query==i, "end_init"]
  chain_1[chain_1$query2==i, "start2"] <- chain_1[chain_1$query2==i, "length2"] - chain_1[chain_1$query2==i, "end_init2"]
  
  chain_1[chain_1$query==i & chain_1$strand=="-", "strand"] <- "temp"
  chain_1[chain_1$query==i & chain_1$strand=="+", "strand"] <- "-"
  chain_1[chain_1$query==i & chain_1$strand=="temp", "strand"] <- "+"
}



#reorient specific chromsomes (y) in specific taxa (bin_id_order[x])

for (i in c(paste(bin_id_order[1], "1", sep = "_"),
            paste(bin_id_order[1], "8", sep = "_"),
            paste(bin_id_order[1], "12", sep = "_"),
            paste(bin_id_order[3], "8", sep = "_"),
            paste(bin_id_order[3], "19", sep = "_"),
            paste(bin_id_order[3], "17", sep = "_"),
            paste(bin_id_order[3], "7", sep = "_"),
            paste(bin_id_order[3], "11", sep = "_"),
            paste(bin_id_order[3], "16", sep = "_"),
            paste(bin_id_order[3], "18", sep = "_"),
            paste(bin_id_order[3], "10", sep = "_"),
            paste(bin_id_order[3], "12", sep = "_"),
            paste(bin_id_order[4], "6", sep = "_"),
            paste(bin_id_order[4], "10", sep = "_"),
            paste(bin_id_order[5], "1", sep = "_"),
            paste(bin_id_order[5], "3", sep = "_"),
            paste(bin_id_order[5], "13", sep = "_"),
            paste(bin_id_order[5], "6", sep = "_"),
            paste(bin_id_order[5], "12", sep = "_"),
            paste(bin_id_order[5], "10", sep = "_"),
            paste(bin_id_order[6], "13", sep = "_"),
            paste(bin_id_order[6], "14", sep = "_"),
            paste(bin_id_order[6], "17", sep = "_"),
            paste(bin_id_order[6], "8", sep = "_"),
            paste(bin_id_order[6], "2", sep = "_"),
            paste(bin_id_order[6], "19", sep = "_"),
            paste(bin_id_order[7], "5", sep = "_"),
            paste(bin_id_order[7], "9", sep = "_"),
            paste(bin_id_order[7], "6", sep = "_"),
            paste(bin_id_order[7], "4", sep = "_"),
            paste(bin_id_order[8], "11", sep = "_"),
            paste(bin_id_order[8], "12", sep = "_"),
            paste(bin_id_order[8], "18", sep = "_"),
            paste(bin_id_order[8], "6", sep = "_"),
            paste(bin_id_order[8], "3", sep = "_"),
            paste(bin_id_order[8], "13", sep = "_"),
            paste(bin_id_order[8], "14", sep = "_"),
            paste(bin_id_order[8], "20", sep = "_"),
            paste(bin_id_order[9], "8", sep = "_"),
            paste(bin_id_order[10], "16", sep = "_"),
            paste(bin_id_order[10], "1", sep = "_"),
            paste(bin_id_order[10], "12", sep = "_"),
            paste(bin_id_order[12], "4", sep = "_"),
            paste(bin_id_order[12], "13", sep = "_"),
            paste(bin_id_order[12], "8", sep = "_"),
            paste(bin_id_order[12], "9", sep = "_"),
            paste(bin_id_order[12], "11", sep = "_")
)) {
    print(i)
   
  chain_1[chain_1$query==i, "end"] <- chain_1[chain_1$query==i, "length"] - chain_1[chain_1$query==i, "start"]
  chain_1[chain_1$query2==i, "end2"] <- chain_1[chain_1$query2==i, "length2"] - chain_1[chain_1$query2==i, "start2"]
  
  chain_1[chain_1$query==i, "start"] <- chain_1[chain_1$query==i, "length"] - chain_1[chain_1$query==i, "end_init"]
  chain_1[chain_1$query2==i, "start2"] <- chain_1[chain_1$query2==i, "length2"] - chain_1[chain_1$query2==i, "end_init2"]
  
    chain_1[chain_1$query==i & chain_1$strand=="-", "strand"] <- "temp"
    chain_1[chain_1$query==i & chain_1$strand=="+", "strand"] <- "-"
    chain_1[chain_1$query==i & chain_1$strand=="temp", "strand"] <- "+"
}

#strand gives highlight in the link in the plot, remove if not wanted
#chain_1$strand <- "+"
link_busco <- chain_1[,c(1:12)]


link_busco$colour_block <- paste("ilMelCinx1", link_busco$colour_block, sep = "_")
write.table(link_busco, file = paste("../intermediate/", prefix, "_link_busco_ordered6.table", sep = ""))

#write.table(sequences, "refined_sequences.table")

########################

link_busco <- read.table(file = paste("../intermediate/", prefix, "_link_busco_ordered6.table", sep = ""))

sequences <- read.table(file = paste("../intermediate/Mec_Mel_sequence_order_mod6.table"))


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

#set colours

synt_col <- rev(as.vector(c(pals::kelly(),pals::polychrome()))[2:c(1+length(unique(link_busco$colour_block)))])
# #change grey
synt_col[9] <-"#85660D"
synt_col[10] <- "#66B0FF"  
synt_col[11] <- "#f4ed00"  
# 
# Make the ribbon plot - these layers can be fully customized as needed!
# 
# If synteny blocks or looking at gene level use geom_link which create polygons instead of lines 
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


synteny_plot <- make_plot(link_busco2, sequences2, add_scale_bar = F)

synteny_plot
# Save the ribbon plot
# 
ggsave(paste("../output/", prefix, "_synteny_gggenomes6", "ms40.png", sep = ""), synteny_plot,
       units = "cm", width = 40, height = 40, bg = "white")

ggsave(paste("../output/", prefix, "_synteny_gggenomes6", "ms40.pdf", sep = ""), synteny_plot,
       units = "cm", width = 40, height = 40, bg = "white")

ggsave(paste("../output/", prefix, "_ynteny_gggenomes6", "ms.png", sep = ""), synteny_plot,
       units = "cm", width = 50, height = 40, bg = "white")

ggsave(paste("../output/", prefix, "_ynteny_gggenomes6", "ms.pdf", sep = ""), synteny_plot,
       units = "cm", width = 50, height = 40, bg = "white")

#ggsave(paste(args$prefix, ".png", sep = ""), synteny_plot,
#       units = "cm", widt = 50, height = 20, bg = "white")




#cat(paste("Plot saved:", paste(args$prefix, ".png", sep = ""), "\n", sep = " "))
