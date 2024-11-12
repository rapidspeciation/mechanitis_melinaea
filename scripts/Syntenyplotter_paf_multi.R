#!/usr/bin/env Rscript --vanilla

#script by Karin Näsvall using syntenyPlotter
#https://github.com/Farre-lab/syntenyPlotteR
#This script uses the output from minimap2 whole genome alignments.
#Creates chain files for syntenyPlotter.
#for multiplotting, do not flip the chromosomes

#install.packages("devtools")
#devtools::install_github("marta-fb/syntenyPlotteR")
library(syntenyPlotteR)
#install.packages("pals")
library(pals)
library(dplyr)
library(ggplot2)

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

paf_table <- read.csv2(table_file, sep = "\t", header = F)[,1:12]

colnames(paf_table) <- c("query", "Qseq_length", "Qstart", "Qend", "strand", "reference", "Rseq_length", "Rstart", "Rend",  "matches", "total","MQ")

head(paf_table)

#check distribution of alignments before filtering 
#hist(paf_table$Rend - paf_table$Rstart, xlim = c(0,100000), breaks = 3000)
#hist(paf_table$total, xlim = c(0,100000), breaks = 3000)

#hist(paf_table[paf_table$total > 100000,"total"])

length(paf_table[paf_table$total > 10000,c("total")])

#filter , used 100-500 kb alignments since closely related taxa, needs to be adjusted
#paf_table <- paf_table[paf_table$total > 500000 & paf_table$MQ==60,]
#paf_table <- paf_table[paf_table$total > 100000 & paf_table$MQ==60,]
paf_table <- paf_table[paf_table$total > 100000,]

#remove small scaffolds
paf_table <- filter(paf_table, !grepl("unloc|SCAFFOLD", query))
paf_table <- filter(paf_table, !grepl("unloc|SCAFFOLD", reference))

paf_table$refID <- taxa1
paf_table$queryID <- taxa2

#correct column order for the chainfile
chain_table <- paf_table[, c("reference", "Rstart", "Rend", "query", "Qstart", "Qend", "strand", "refID", "queryID", "Rseq_length", "Qseq_length")]


#change names of chr OBS customise!!!
unique(chain_table$reference)
unique(chain_table$query)
#get default colours if chr names are just numbers
#otherwise have to assign colours manually for the plot
chain_table$reference <- sub("SUPER_", "", chain_table$reference)
chain_table$query <- sub("SUPER_", "", chain_table$query)

chain_table$reference <- sub(paste("il", taxa1, "_", sep = ""), "", chain_table$reference)
chain_table$query <- sub(paste("il", taxa1, "_", sep = ""), "", chain_table$query)
head(chain_table)

#order df after chr length in ref taxa, will also order the other taxa after first appearance in the chain file
chain_table <- 
  chain_table %>% arrange(desc(Rseq_length))


#make size_file for input (and order of chr in the plot)

chr_size_q <- unique(chain_table[,c("query", "Qseq_length", "queryID")])
colnames(chr_size_q) <- c("chr", "seq_length", "ID")
chr_size_r <- unique(chain_table[,c("reference", "Rseq_length", "refID")])
colnames(chr_size_r) <- c("chr", "seq_length", "ID")
chr_size <- rbind(chr_size_q, chr_size_r)


#writing input files
names(chr_size) <- NULL
write.table(chr_size, file = paste("intermediate/chr_length_", taxa2, taxa1, ".txt", sep = ""), sep = "\t", row.names = F)

names(chain_table) <- NULL
write.table(chain_table, file = paste("intermediate/chain_", taxa1, taxa2, ".txt", sep = ""), sep = "\t", row.names = F)

#produce the synteny figures
draw.linear(directory = "plots", output = paste("synt_", taxa1, taxa2, Sys.Date(), sep = ""), 
            paste("intermediate/chr_length_", taxa2, taxa1, ".txt", sep = ""), 
            paste("intermediate/chain_", taxa1, taxa2, ".txt", sep = ""),  
            fileformat = "pdf", w=13,h=5)

# ggsave(filename = paste("plots/synt_", taxa1, taxa2, Sys.Date(), ".pdf", sep = ""), 
#        device = "png", width = 13, height = 5)

