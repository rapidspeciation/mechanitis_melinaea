#!/usr/bin/env Rscript --vanilla

#script by Karin Näsvall using syntenyPlotter
#https://github.com/Farre-lab/syntenyPlotteR
#This script uses the output from minimap2 whole genome alignments.
#Creates chain files for syntenyPlotter.

#install.packages("devtools")
#devtools::install_github("marta-fb/syntenyPlotteR")
library(syntenyPlotteR)
#install.packages("pals")
library(pals)
library(dplyr)
library(ggplot2)
library(tidyverse)
#rm(list = ls())
#
#dir organisation: intermediate/, plots/


######
####If using command line and argument uncomment these (and comment away the variable lines below):
cmd_args <- commandArgs(trailingOnly = TRUE)
table_file <- cmd_args[1]
taxa1=cmd_args[2]
taxa2=cmd_args[3]

#################
# Variables if using hardcoding, change these variables:
#output from minimap2 paf-format
# table_file <- ""
# #reference (the taxa in the sixth column of the output in minimap2)
# taxa1=""
# #query (the taxa in the first column of the output in minimap2)
# taxa2=""

#################

# list_chr_files <- read.table("output_species_comp/list_chr_files.txt")
# 
# list_chr_files$V1 <- paste("output_species_comp/", list_chr_files$V1, sep="")
# 
# for(file in list_chr_files$V1){
#   chr_size <- read_delim(list_chr_files$V1, col_names = FALSE)
# }
# 
# 
# 
# chr_size <- read_delim(list_chr_files$V1, col_names = FALSE)
# 
# chr_size$X3 <- as.factor(chr_size$X3)
# split.df <- split(chr_size, f = chr_size$X3, drop = TRUE)
# 

#################
#################
#flipping chr

chain_table <- read.table("output_species_comp/intermediate/chain_MecLysi1MecMaza1.txt")[,1:11]
colnames(chain_table) <- c("reference", "Rstart", "Rend", "query", "Qstart", "Qend", "strand", "refID", "queryID", "Rseq_length", "Qseq_length")

taxa1 <- unique(chain_table$refID)
taxa2 <- unique(chain_table$queryID)

chain_1 <- chain_table
chain_1$Rstart_init <- chain_1$Rstart
chain_1$Rend_init <- chain_1$Rend

chain_1$Qstart_init <- chain_1$Qstart
chain_1$Qend_init <- chain_1$Qend


for (i in c("1", "W1", "Z1")) {
    
    chain_1[chain_1$query==i, "Qend"] <- chain_1[chain_1$query==i, "Qseq_length"] - chain_1[chain_1$query==i, "Qstart"]
    
    chain_1[chain_1$query==i, "Qstart"] <- chain_1[chain_1$query==i, "Qseq_length"] - chain_1[chain_1$query==i, "Qend_init"]
    
    chain_1[chain_1$query==i & chain_1$strand=="-", "strand"] <- "temp"
    chain_1[chain_1$query==i & chain_1$strand=="+", "strand"] <- "-"
    chain_1[chain_1$query==i & chain_1$strand=="temp", "strand"] <- "+"
}

for (i in c("10", "14", "11","21","2","16", "Z1")) {
  
  chain_1[chain_1$reference==i, "Rend"] <- chain_1[chain_1$reference==i, "Rseq_length"] - chain_1[chain_1$reference==i, "Rstart"]
  
  chain_1[chain_1$reference==i, "Rstart"] <- chain_1[chain_1$reference==i, "Rseq_length"] - chain_1[chain_1$reference==i, "Rend_init"]
  
  chain_1[chain_1$reference==i & chain_1$strand=="-", "strand"] <- "temp"
  chain_1[chain_1$reference==i & chain_1$strand=="+", "strand"] <- "-"
  chain_1[chain_1$reference==i & chain_1$strand=="temp", "strand"] <- "+"
}



chain_table <- chain_1


names(chain_table) <- NULL
write.table(chain_table, file = paste("output_species_comp/intermediate/chain_", taxa1, taxa2, "_reordered.txt", sep = ""), sep = "\t", row.names = F)

#####
#flipping chr

chain_table <- read.table("output_species_comp/intermediate/chain_MecMaza1MecMena1.txt")[,1:11]
colnames(chain_table) <- c("reference", "Rstart", "Rend", "query", "Qstart", "Qend", "strand", "refID", "queryID", "Rseq_length", "Qseq_length")

taxa1 <- unique(chain_table$refID)
taxa2 <- unique(chain_table$queryID)

chain_1 <- chain_table
chain_1$Rstart_init <- chain_1$Rstart
chain_1$Rend_init <- chain_1$Rend

chain_1$Qstart_init <- chain_1$Qstart
chain_1$Qend_init <- chain_1$Qend


for (i in c("6", "W","Z", "11")) {
  
  chain_1[chain_1$query==i, "Qend"] <- chain_1[chain_1$query==i, "Qseq_length"] - chain_1[chain_1$query==i, "Qstart"]
  
  chain_1[chain_1$query==i, "Qstart"] <- chain_1[chain_1$query==i, "Qseq_length"] - chain_1[chain_1$query==i, "Qend_init"]
  
  chain_1[chain_1$query==i & chain_1$strand=="-", "strand"] <- "temp"
  chain_1[chain_1$query==i & chain_1$strand=="+", "strand"] <- "-"
  chain_1[chain_1$query==i & chain_1$strand=="temp", "strand"] <- "+"
}


for (i in c("1", "Z1")) {
  
  chain_1[chain_1$reference==i, "Rend"] <- chain_1[chain_1$reference==i, "Rseq_length"] - chain_1[chain_1$reference==i, "Rstart"]
  
  chain_1[chain_1$reference==i, "Rstart"] <- chain_1[chain_1$reference==i, "Rseq_length"] - chain_1[chain_1$reference==i, "Rend_init"]
  
  chain_1[chain_1$reference==i & chain_1$strand=="-", "strand"] <- "temp"
  chain_1[chain_1$reference==i & chain_1$strand=="+", "strand"] <- "-"
  chain_1[chain_1$reference==i & chain_1$strand=="temp", "strand"] <- "+"
}



chain_table <- chain_1


names(chain_table) <- NULL
write.table(chain_table, file = paste("output_species_comp/intermediate/chain_", taxa1, taxa2, "_reordered.txt", sep = ""), sep = "\t", row.names = F)

#####

#####
#flipping chr

chain_table <- read.table("output_species_comp/intermediate/chain_MecMena1MecMess1.txt")[,1:11]
colnames(chain_table) <- c("reference", "Rstart", "Rend", "query", "Qstart", "Qend", "strand", "refID", "queryID", "Rseq_length", "Qseq_length")

taxa1 <- unique(chain_table$refID)
taxa2 <- unique(chain_table$queryID)

chain_1 <- chain_table
chain_1$Rstart_init <- chain_1$Rstart
chain_1$Rend_init <- chain_1$Rend

chain_1$Qstart_init <- chain_1$Qstart
chain_1$Qend_init <- chain_1$Qend


for (i in c("W", "Z")) {
  
  chain_1[chain_1$query==i, "Qend"] <- chain_1[chain_1$query==i, "Qseq_length"] - chain_1[chain_1$query==i, "Qstart"]
  
  chain_1[chain_1$query==i, "Qstart"] <- chain_1[chain_1$query==i, "Qseq_length"] - chain_1[chain_1$query==i, "Qend_init"]
  
  chain_1[chain_1$query==i & chain_1$strand=="-", "strand"] <- "temp"
  chain_1[chain_1$query==i & chain_1$strand=="+", "strand"] <- "-"
  chain_1[chain_1$query==i & chain_1$strand=="temp", "strand"] <- "+"
}


for (i in c("6", "W", "Z", "11")) {
  
  chain_1[chain_1$reference==i, "Rend"] <- chain_1[chain_1$reference==i, "Rseq_length"] - chain_1[chain_1$reference==i, "Rstart"]
  
  chain_1[chain_1$reference==i, "Rstart"] <- chain_1[chain_1$reference==i, "Rseq_length"] - chain_1[chain_1$reference==i, "Rend_init"]
  
  chain_1[chain_1$reference==i & chain_1$strand=="-", "strand"] <- "temp"
  chain_1[chain_1$reference==i & chain_1$strand=="+", "strand"] <- "-"
  chain_1[chain_1$reference==i & chain_1$strand=="temp", "strand"] <- "+"
}



chain_table <- chain_1


names(chain_table) <- NULL
write.table(chain_table, file = paste("output_species_comp/intermediate/chain_", taxa1, taxa2, "_reordered.txt", sep = ""), sep = "\t", row.names = F)

#####


#####
#flipping chr

chain_table <- read.table("output_species_comp/intermediate/chain_MecMess1MecPoly1.txt")[,1:11]
colnames(chain_table) <- c("reference", "Rstart", "Rend", "query", "Qstart", "Qend", "strand", "refID", "queryID", "Rseq_length", "Qseq_length")

taxa1 <- unique(chain_table$refID)
taxa2 <- unique(chain_table$queryID)

chain_1 <- chain_table
chain_1$Rstart_init <- chain_1$Rstart
chain_1$Rend_init <- chain_1$Rend

chain_1$Qstart_init <- chain_1$Qstart
chain_1$Qend_init <- chain_1$Qend


for (i in c("14", "W", "13", "9", "Z")) {
  
  chain_1[chain_1$query==i, "Qend"] <- chain_1[chain_1$query==i, "Qseq_length"] - chain_1[chain_1$query==i, "Qstart"]
  
  chain_1[chain_1$query==i, "Qstart"] <- chain_1[chain_1$query==i, "Qseq_length"] - chain_1[chain_1$query==i, "Qend_init"]
  
  chain_1[chain_1$query==i & chain_1$strand=="-", "strand"] <- "temp"
  chain_1[chain_1$query==i & chain_1$strand=="+", "strand"] <- "-"
  chain_1[chain_1$query==i & chain_1$strand=="temp", "strand"] <- "+"
}


for (i in c("W", "Z")) {
  
  chain_1[chain_1$reference==i, "Rend"] <- chain_1[chain_1$reference==i, "Rseq_length"] - chain_1[chain_1$reference==i, "Rstart"]
  
  chain_1[chain_1$reference==i, "Rstart"] <- chain_1[chain_1$reference==i, "Rseq_length"] - chain_1[chain_1$reference==i, "Rend_init"]
  
  chain_1[chain_1$reference==i & chain_1$strand=="-", "strand"] <- "temp"
  chain_1[chain_1$reference==i & chain_1$strand=="+", "strand"] <- "-"
  chain_1[chain_1$reference==i & chain_1$strand=="temp", "strand"] <- "+"
}



chain_table <- chain_1


names(chain_table) <- NULL
write.table(chain_table, file = paste("output_species_comp/intermediate/chain_", taxa1, taxa2, "_reordered.txt", sep = ""), sep = "\t", row.names = F)

#####


#####
#flipping chr

chain_table <- read.table("output_species_comp/intermediate/chain_MelLudo1MelIsoc1.txt")[,1:11]
colnames(chain_table) <- c("reference", "Rstart", "Rend", "query", "Qstart", "Qend", "strand", "refID", "queryID", "Rseq_length", "Qseq_length")

taxa1 <- unique(chain_table$refID)
taxa2 <- unique(chain_table$queryID)

#change 2 to Z in isocomma
chain_table[chain_table$query==2,"query"] <- "Z"


chain_1 <- chain_table
chain_1$Rstart_init <- chain_1$Rstart
chain_1$Rend_init <- chain_1$Rend

chain_1$Qstart_init <- chain_1$Qstart
chain_1$Qend_init <- chain_1$Qend

for (i in c("14","12","6")) {
  
  chain_1[chain_1$query==i, "Qend"] <- chain_1[chain_1$query==i, "Qseq_length"] - chain_1[chain_1$query==i, "Qstart"]
  
  chain_1[chain_1$query==i, "Qstart"] <- chain_1[chain_1$query==i, "Qseq_length"] - chain_1[chain_1$query==i, "Qend_init"]
  
  chain_1[chain_1$query==i & chain_1$strand=="-", "strand"] <- "temp"
  chain_1[chain_1$query==i & chain_1$strand=="+", "strand"] <- "-"
  chain_1[chain_1$query==i & chain_1$strand=="temp", "strand"] <- "+"
}


for (i in c("14", "4","19", "5","9","15","6")) {
  
  chain_1[chain_1$reference==i, "Rend"] <- chain_1[chain_1$reference==i, "Rseq_length"] - chain_1[chain_1$reference==i, "Rstart"]
  
  chain_1[chain_1$reference==i, "Rstart"] <- chain_1[chain_1$reference==i, "Rseq_length"] - chain_1[chain_1$reference==i, "Rend_init"]
  
  chain_1[chain_1$reference==i & chain_1$strand=="-", "strand"] <- "temp"
  chain_1[chain_1$reference==i & chain_1$strand=="+", "strand"] <- "-"
  chain_1[chain_1$reference==i & chain_1$strand=="temp", "strand"] <- "+"
}



chain_table <- chain_1


names(chain_table) <- NULL
write.table(chain_table, file = paste("output_species_comp/intermediate/chain_", taxa1, taxa2, "_reordered.txt", sep = ""), sep = "\t", row.names = F)

#####

#####
#flipping chr

chain_table <- read.table("output_species_comp/intermediate/chain_MelIsoc1MelMoth8.txt")[,1:11]
colnames(chain_table) <- c("reference", "Rstart", "Rend", "query", "Qstart", "Qend", "strand", "refID", "queryID", "Rseq_length", "Qseq_length")

taxa1 <- unique(chain_table$refID)
taxa2 <- unique(chain_table$queryID)

#change 2 to Z in isocomma
chain_table[chain_table$reference==2,"reference"] <- "Z"


chain_1 <- chain_table
chain_1$Rstart_init <- chain_1$Rstart
chain_1$Rend_init <- chain_1$Rend

chain_1$Qstart_init <- chain_1$Qstart
chain_1$Qend_init <- chain_1$Qend


for (i in c("5","6","W1", "Z2", "7", "Z1", "9")) {
  
  chain_1[chain_1$query==i, "Qend"] <- chain_1[chain_1$query==i, "Qseq_length"] - chain_1[chain_1$query==i, "Qstart"]
  
  chain_1[chain_1$query==i, "Qstart"] <- chain_1[chain_1$query==i, "Qseq_length"] - chain_1[chain_1$query==i, "Qend_init"]
  
  chain_1[chain_1$query==i & chain_1$strand=="-", "strand"] <- "temp"
  chain_1[chain_1$query==i & chain_1$strand=="+", "strand"] <- "-"
  chain_1[chain_1$query==i & chain_1$strand=="temp", "strand"] <- "+"
}

for (i in c("14","12","9")) {
  
  chain_1[chain_1$reference==i, "Rend"] <- chain_1[chain_1$reference==i, "Rseq_length"] - chain_1[chain_1$reference==i, "Rstart"]
  
  chain_1[chain_1$reference==i, "Rstart"] <- chain_1[chain_1$reference==i, "Rseq_length"] - chain_1[chain_1$reference==i, "Rend_init"]
  
  chain_1[chain_1$reference==i & chain_1$strand=="-", "strand"] <- "temp"
  chain_1[chain_1$reference==i & chain_1$strand=="+", "strand"] <- "-"
  chain_1[chain_1$reference==i & chain_1$strand=="temp", "strand"] <- "+"
}



chain_table <- chain_1


names(chain_table) <- NULL
write.table(chain_table, file = paste("output_species_comp/intermediate/chain_", taxa1, taxa2, "_reordered.txt", sep = ""), sep = "\t", row.names = F)

#####

#####
#flipping chr

chain_table <- read.table("output_species_comp/intermediate/chain_MelMoth8MelMeno1.txt")[,1:11]
colnames(chain_table) <- c("reference", "Rstart", "Rend", "query", "Qstart", "Qend", "strand", "refID", "queryID", "Rseq_length", "Qseq_length")

taxa1 <- unique(chain_table$refID)
taxa2 <- unique(chain_table$queryID)


chain_1 <- chain_table
chain_1$Rstart_init <- chain_1$Rstart
chain_1$Rend_init <- chain_1$Rend

chain_1$Qstart_init <- chain_1$Qstart
chain_1$Qend_init <- chain_1$Qend


for (i in c("8", "7","15", "1", "11", "6", "16","2")) {
  
  chain_1[chain_1$query==i, "Qend"] <- chain_1[chain_1$query==i, "Qseq_length"] - chain_1[chain_1$query==i, "Qstart"]
  
  chain_1[chain_1$query==i, "Qstart"] <- chain_1[chain_1$query==i, "Qseq_length"] - chain_1[chain_1$query==i, "Qend_init"]
  
  chain_1[chain_1$query==i & chain_1$strand=="-", "strand"] <- "temp"
  chain_1[chain_1$query==i & chain_1$strand=="+", "strand"] <- "-"
  chain_1[chain_1$query==i & chain_1$strand=="temp", "strand"] <- "+"
}

for (i in c("10", "9", "4", "Z1", "W1")) {
  
  chain_1[chain_1$reference==i, "Rend"] <- chain_1[chain_1$reference==i, "Rseq_length"] - chain_1[chain_1$reference==i, "Rstart"]
  
  chain_1[chain_1$reference==i, "Rstart"] <- chain_1[chain_1$reference==i, "Rseq_length"] - chain_1[chain_1$reference==i, "Rend_init"]
  
  chain_1[chain_1$reference==i & chain_1$strand=="-", "strand"] <- "temp"
  chain_1[chain_1$reference==i & chain_1$strand=="+", "strand"] <- "-"
  chain_1[chain_1$reference==i & chain_1$strand=="temp", "strand"] <- "+"
}



chain_table <- chain_1


names(chain_table) <- NULL
write.table(chain_table, file = paste("output_species_comp/intermediate/chain_", taxa1, taxa2, "_reordered.txt", sep = ""), sep = "\t", row.names = F)

#####

#####
#flipping chr

chain_table <- read.table("output_species_comp/intermediate/chain_MelMeno1MelMars1.txt")[,1:11]
colnames(chain_table) <- c("reference", "Rstart", "Rend", "query", "Qstart", "Qend", "strand", "refID", "queryID", "Rseq_length", "Qseq_length")

taxa1 <- unique(chain_table$refID)
taxa2 <- unique(chain_table$queryID)


chain_1 <- chain_table
chain_1$Rstart_init <- chain_1$Rstart
chain_1$Rend_init <- chain_1$Rend

chain_1$Qstart_init <- chain_1$Qstart
chain_1$Qend_init <- chain_1$Qend


for (i in c("6","7", "8", "Z" )) {
  
  chain_1[chain_1$query==i, "Qend"] <- chain_1[chain_1$query==i, "Qseq_length"] - chain_1[chain_1$query==i, "Qstart"]
  
  chain_1[chain_1$query==i, "Qstart"] <- chain_1[chain_1$query==i, "Qseq_length"] - chain_1[chain_1$query==i, "Qend_init"]
  
  chain_1[chain_1$query==i & chain_1$strand=="-", "strand"] <- "temp"
  chain_1[chain_1$query==i & chain_1$strand=="+", "strand"] <- "-"
  chain_1[chain_1$query==i & chain_1$strand=="temp", "strand"] <- "+"
}

for (i in c( "14", "7","15","4", "11", "18", "5", "16", "19", "12")) {
  
  chain_1[chain_1$reference==i, "Rend"] <- chain_1[chain_1$reference==i, "Rseq_length"] - chain_1[chain_1$reference==i, "Rstart"]
  
  chain_1[chain_1$reference==i, "Rstart"] <- chain_1[chain_1$reference==i, "Rseq_length"] - chain_1[chain_1$reference==i, "Rend_init"]
  
  chain_1[chain_1$reference==i & chain_1$strand=="-", "strand"] <- "temp"
  chain_1[chain_1$reference==i & chain_1$strand=="+", "strand"] <- "-"
  chain_1[chain_1$reference==i & chain_1$strand=="temp", "strand"] <- "+"
}



chain_table <- chain_1


names(chain_table) <- NULL
write.table(chain_table, file = paste("output_species_comp/intermediate/chain_", taxa1, taxa2, "_reordered.txt", sep = ""), sep = "\t", row.names = F)

#####

#write.table(chr_size, file = "output_species_comp/intermediate/chr_length_all.txt", col.names = FALSE)
#produce the synteny figures
draw.linear(directory = "output_species_comp/plots", output = paste("synt_all", Sys.Date(),format(Sys.time(), "%X"), sep = "_"), 
            sizefile = "output_species_comp/intermediate/all_chr_length_all_ordered.txt",
            "output_species_comp/intermediate/chain_MecMess1MecPoly1.txt",
            "output_species_comp/intermediate/chain_MecMena1MecMess1.txt",
            "output_species_comp/intermediate/chain_MecMaza1MecMena1.txt",
            "output_species_comp/intermediate/chain_MecLysi1MecMaza1.txt",
            "output_species_comp/intermediate/chain_MelMars1MecLysi1.txt",
            "output_species_comp/intermediate/chain_MelMeno1MelMars1.txt",
            "output_species_comp/intermediate/chain_MelMoth8MelMeno1.txt",
            "output_species_comp/intermediate/chain_MelIsoc1MelMoth8.txt",
            "output_species_comp/intermediate/chain_MelLudo1MelIsoc1.txt",
            fileformat = "pdf", w=13,h=13, colours = c("Z"= "black", "Z1"="black", "Z2"="black", "W"="black", "W1" = "black", "W2"="black", "W3" = "black" ))


#write.table(chr_size, file = "output_species_comp/intermediate/chr_length_all.txt", col.names = FALSE)
#produce the synteny figures
draw.linear(directory = "output_species_comp/plots", output = paste("synt_all", Sys.Date(),format(Sys.time(), "%X"),"_darkgreen", sep = "_"), 
            sizefile = "output_species_comp/intermediate/all_chr_length_all_ordered.txt",
            "output_species_comp/intermediate/chain_MecMess1MecPoly1_reordered.txt",
            "output_species_comp/intermediate/chain_MecMena1MecMess1_reordered.txt",
            "output_species_comp/intermediate/chain_MecMaza1MecMena1_reordered.txt",
            "output_species_comp/intermediate/chain_MecLysi1MecMaza1_reordered.txt",
            "output_species_comp/intermediate/chain_MelMars1MecLysi1.txt",
            "output_species_comp/intermediate/chain_MelMeno1MelMars1_reordered.txt",
            "output_species_comp/intermediate/chain_MelMoth8MelMeno1_reordered.txt",
            "output_species_comp/intermediate/chain_MelIsoc1MelMoth8_reordered.txt",
            "output_species_comp/intermediate/chain_MelLudo1MelIsoc1_reordered.txt",
 fileformat = "pdf", w=7.5,h=6, colours = c("1" = "darkgreen",
                                          "10" = "darkgreen",
                                          "11" = "darkgreen",
                                          "12" = "darkgreen",
                                          "13" = "darkgreen",
                                          "14" = "darkgreen",
                                          "15" = "darkgreen",
                                          "16" = "darkgreen",
                                          "17" = "darkgreen",
                                          "18" = "darkgreen",
                                          "19" = "darkgreen",
                                          "2" = "darkgreen",
                                          "20" = "darkgreen",
                                          "21" = "darkgreen",
                                          "3" = "darkgreen",
                                          "4" = "darkgreen",
                                          "5" = "darkgreen",
                                          "6" = "darkgreen",
                                          "7" = "darkgreen",
                                          "8" = "darkgreen",
                                          "9" = "darkgreen",
                                          "W" = "darkgreen",
                                          "W1" = "darkgreen",
                                          "W2" = "darkgreen",
                                          "W3" = "darkgreen",
                                          "Z" = "darkgreen",
                                          "Z1" = "darkgreen",
                                          "Z2" = "darkgreen"
                              ))


draw.linear(directory = "output_species_comp/plots", output = paste("synt_all", Sys.Date(),format(Sys.time(), "%X"), "_golden", sep = "_"), 
            sizefile = "output_species_comp/intermediate/all_chr_length_all_ordered.txt",
            "output_species_comp/intermediate/chain_MecMess1MecPoly1_reordered.txt",
            "output_species_comp/intermediate/chain_MecMena1MecMess1_reordered.txt",
            "output_species_comp/intermediate/chain_MecMaza1MecMena1_reordered.txt",
            "output_species_comp/intermediate/chain_MecLysi1MecMaza1_reordered.txt",
            "output_species_comp/intermediate/chain_MelMars1MecLysi1.txt",
            "output_species_comp/intermediate/chain_MelMeno1MelMars1_reordered.txt",
            "output_species_comp/intermediate/chain_MelMoth8MelMeno1_reordered.txt",
            "output_species_comp/intermediate/chain_MelIsoc1MelMoth8_reordered.txt",
            "output_species_comp/intermediate/chain_MelLudo1MelIsoc1_reordered.txt",
            fileformat = "pdf", w=12.5,h=10, opacity = 0.7, colours = c("1" = "goldenrod3",
                                                     "10" = "goldenrod3",
                                                     "11" = "goldenrod3",
                                                     "12" = "goldenrod3",
                                                     "13" = "goldenrod3",
                                                     "14" = "goldenrod3",
                                                     "15" = "goldenrod3",
                                                     "16" = "goldenrod3",
                                                     "17" = "goldenrod3",
                                                     "18" = "goldenrod3",
                                                     "19" = "goldenrod3",
                                                     "2" = "goldenrod3",
                                                     "20" = "goldenrod3",
                                                     "21" = "goldenrod3",
                                                     "3" = "goldenrod3",
                                                     "4" = "goldenrod3",
                                                     "5" = "goldenrod3",
                                                     "6" = "goldenrod3",
                                                     "7" = "goldenrod3",
                                                     "8" = "goldenrod3",
                                                     "9" = "goldenrod3",
                                                     "W" = "goldenrod3",
                                                     "W1" = "goldenrod3",
                                                     "W2" = "goldenrod3",
                                                     "W3" = "goldenrod3",
                                                     "Z" = "goldenrod3",
                                                     "Z1" = "goldenrod3",
                                                     "Z2" = "goldenrod3"
            ))


draw.linear(directory = "output_species_comp/plots", output = paste("synt_all", Sys.Date(),format(Sys.time(), "%X"),"_royalblue2", sep = "_"), 
            sizefile = "output_species_comp/intermediate/all_chr_length_all_ordered.txt",
            "output_species_comp/intermediate/chain_MecMess1MecPoly1_reordered.txt",
            "output_species_comp/intermediate/chain_MecMena1MecMess1_reordered.txt",
            "output_species_comp/intermediate/chain_MecMaza1MecMena1_reordered.txt",
            "output_species_comp/intermediate/chain_MecLysi1MecMaza1_reordered.txt",
            "output_species_comp/intermediate/chain_MelMars1MecLysi1.txt",
            "output_species_comp/intermediate/chain_MelMeno1MelMars1_reordered.txt",
            "output_species_comp/intermediate/chain_MelMoth8MelMeno1_reordered.txt",
            "output_species_comp/intermediate/chain_MelIsoc1MelMoth8_reordered.txt",
            "output_species_comp/intermediate/chain_MelLudo1MelIsoc1_reordered.txt",
            fileformat = "pdf", w=8,h=6, colours = c("1" = "royalblue2",
                                                      "10" = "royalblue2",
                                                      "11" = "royalblue2",
                                                      "12" = "royalblue2",
                                                      "13" = "royalblue2",
                                                      "14" = "royalblue2",
                                                      "15" = "royalblue2",
                                                      "16" = "royalblue2",
                                                      "17" = "royalblue2",
                                                      "18" = "royalblue2",
                                                      "19" = "royalblue2",
                                                      "2" = "royalblue2",
                                                      "20" = "royalblue2",
                                                      "21" = "royalblue2",
                                                      "3" = "royalblue2",
                                                      "4" = "royalblue2",
                                                      "5" = "royalblue2",
                                                      "6" = "royalblue2",
                                                      "7" = "royalblue2",
                                                      "8" = "royalblue2",
                                                      "9" = "royalblue2",
                                                      "W" = "royalblue2",
                                                      "W1" = "royalblue2",
                                                      "W2" = "royalblue2",
                                                      "W3" = "royalblue2",
                                                      "Z" = "royalblue2",
                                                      "Z1" = "royalblue2",
                                                      "Z2" = "royalblue2"
            ))


draw.linear(directory = "output_species_comp/plots", output = paste("synt_all", Sys.Date(),format(Sys.time(), "%X"), "_golden", sep = "_"), 
            sizefile = "output_species_comp/intermediate/all_chr_length_all_ordered.txt",
            "output_species_comp/intermediate/chain_MecMess1MecPoly1_reordered.txt",
            "output_species_comp/intermediate/chain_MecMena1MecMess1_reordered.txt",
            "output_species_comp/intermediate/chain_MecMaza1MecMena1_reordered.txt",
            "output_species_comp/intermediate/chain_MecLysi1MecMaza1_reordered.txt",
            "output_species_comp/intermediate/chain_MelMars1MecLysi1.txt",
            "output_species_comp/intermediate/chain_MelMeno1MelMars1_reordered.txt",
            "output_species_comp/intermediate/chain_MelMoth8MelMeno1_reordered.txt",
            "output_species_comp/intermediate/chain_MelIsoc1MelMoth8_reordered.txt",
            "output_species_comp/intermediate/chain_MelLudo1MelIsoc1_reordered.txt",
            fileformat = "pdf", w=7.5,h=6, opacity = 0.5, colours = c("1" = "goldenrod3",
                                                                      "10" = "goldenrod3",
                                                                      "11" = "goldenrod3",
                                                                      "12" = "goldenrod3",
                                                                      "13" = "goldenrod3",
                                                                      "14" = "goldenrod3",
                                                                      "15" = "goldenrod3",
                                                                      "16" = "goldenrod3",
                                                                      "17" = "goldenrod3",
                                                                      "18" = "goldenrod3",
                                                                      "19" = "goldenrod3",
                                                                      "2" = "goldenrod3",
                                                                      "20" = "goldenrod3",
                                                                      "21" = "goldenrod3",
                                                                      "3" = "goldenrod3",
                                                                      "4" = "goldenrod3",
                                                                      "5" = "goldenrod3",
                                                                      "6" = "goldenrod3",
                                                                      "7" = "goldenrod3",
                                                                      "8" = "goldenrod3",
                                                                      "9" = "goldenrod3",
                                                                      "W" = "goldenrod3",
                                                                      "W1" = "goldenrod3",
                                                                      "W2" = "goldenrod3",
                                                                      "W3" = "goldenrod3",
                                                                      "Z" = "goldenrod3",
                                                                      "Z1" = "goldenrod3",
                                                                      "Z2" = "goldenrod3"
            ))


# ggsave(filename = paste("plots/synt_", taxa1, taxa2, Sys.Date(), ".pdf", sep = ""), 
#        device = "png", width = 13, height = 5)
#c("Z"= "black", "Z1"="black", "Z2"="black", "W"="black", "W1" = "black", "W2"="black", "W3" = "black" )
#colours = c("1" = "yellow3", "2" = "yellow3", "14" = "yellow3", "16" = "yellow3", "18" = "yellow3", "7" = "yellow3", "10" = "yellow3")
#
#
draw.linear(directory = "output_species_comp/plots", output = paste("synt_mec", Sys.Date(),format(Sys.time(), "%X"), sep = "_"), 
            sizefile = "output_species_comp/intermediate/all_chr_length_all_ordered.txt",
            "output_species_comp/intermediate/chain_MecMess1MecPoly1_reordered.txt",
            "output_species_comp/intermediate/chain_MecMena1MecMess1_reordered.txt",
            "output_species_comp/intermediate/chain_MecMaza1MecMena1_reordered.txt",
            "output_species_comp/intermediate/chain_MecLysi1MecMaza1_reordered.txt",
            fileformat = "pdf", w=10,h=5, colours = c("49"="grey"))

draw.linear(directory = "output_species_comp/plots", output = paste("synt_mel", Sys.Date(),format(Sys.time(), "%X"), sep = "_"), 
            sizefile = "output_species_comp/intermediate/all_chr_length_all_ordered.txt",
            "output_species_comp/intermediate/chain_MelMeno1MelMars1_reordered.txt",
            "output_species_comp/intermediate/chain_MelMoth8MelMeno1_reordered.txt",
            "output_species_comp/intermediate/chain_MelIsoc1MelMoth8_reordered.txt",
            "output_species_comp/intermediate/chain_MelLudo1MelIsoc1_reordered.txt",
            fileformat = "pdf", w=10,h=5, colours = c("49"="grey"))

