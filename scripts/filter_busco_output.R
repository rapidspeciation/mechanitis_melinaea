#this script renames the taxa and chromosomes, and filter so only chromosome scaffolds remain from the output of busco full_table.tsv

######################
##Set variables
prefix <- "Mec_Mel"

#busco output full_table.tsv, concatenated for all species of interest, header lines removed
input_link <- "../input/full_table_all_mod.tsv"

#concatenated fai (from samtools faidx) 
input_chr_length <- "../input/chr_length_all.tsv"


######################
#create the sequence file, to get the length of the chromosomes
#read in length of chromosomes
chr_length <- read.csv(input_chr_length, sep = "\t", header = F)
colnames(chr_length) <- c("seq_id", "length", "bin_id")

#make taxa (bin) and adjust names to match busco
chr_length$bin_id <-  sub(".*/", "", chr_length$bin_id)
chr_length$bin_id <-  sub("\\.fai", "", chr_length$bin_id)
chr_length$bin_id <-  sub("\\.gz", "", chr_length$bin_id)
unique(chr_length$bin_id)


######################
#create the links
#links_busco_prel <- read.csv(args$links,
#                             sep = "\t", header = TRUE)

#read in full_table from BUSCO
links_busco_prel <- read.csv(input_link, sep = "\t", header = F)

#only primary
unique(links_busco_prel$V7)

#losing lots of buscos on Z2, duplicated on W
add_Z2 <- links_busco_prel[links_busco_prel$V2=="Duplicated" & links_busco_prel$V3=="SUPER_Z2", c(1,3:7)]

#remove those that are duplicated on Z2
add_Z2 <- 
  add_Z2 %>%
  group_by(V1) %>%
  filter(n()==1) %>%
  ungroup()


#filter only Complete genes in the other taxa
links_busco_prel <- links_busco_prel[links_busco_prel$V2=="Complete", c(1,3:7)]

#add dupl on Z2
links_busco_prel <- rbind(links_busco_prel,add_Z2)

links_busco_prel$V7 <- sub("busco/", "", links_busco_prel$V7)

colnames(links_busco_prel) <- c("block_id","seq_id","start","end", "strand", "bin_id")


links_busco_prel %>%
  group_by(bin_id) %>%
  count()

#add chr length
links_busco_prel <- 
  left_join(links_busco_prel, chr_length)

#filter to use only primary assembly (remove HAP2)
links_busco_prel <- 
  links_busco_prel %>%
  filter(!(str_detect(bin_id, "HAP2") | str_detect(bin_id, "alternative"))) %>%
  filter(!bin_id=="ilMelIsoc1.2.primary.curated.fa")


#refine taxa name
unique(links_busco_prel$bin_id)
links_busco_prel[links_busco_prel$bin_id=="GCF_009731565.1_Dplex_v4_genomic.fa", "bin_id"]="ilDanPlex4"
links_busco_prel[links_busco_prel$bin_id=="GCF_905220565.1_ilMelCinx1.1_genomic.fa", "bin_id"]="ilMelCinx1"
links_busco_prel$bin_id <- sub("\\..*", "", links_busco_prel$bin_id)
links_busco_prel$bin_id <- sub("_.*", "", links_busco_prel$bin_id)

#filter sequences, only keep chromosome level scaffold
links_busco_prel <- links_busco_prel %>%
  filter(!str_detect(seq_id, "SCAFFOLD")) %>%
  filter(!str_detect(seq_id, "unloc")) %>%
  filter(!str_detect(seq_id, "CAUF")) %>%
  filter(!str_detect(seq_id, "MT")) %>%
  filter(!str_detect(seq_id, "NW"))

unique(links_busco_prel$seq_id)

#for assemblies with repository number of the sequences 
#grep ">" input_fa/data_mec/ilMecMess1.1.fasta |awk '{print $1, $7}' |sed 's/>//' >> chr_name_change.list
chr_change <- read.table("../input/chr_name_change.list")


for (i in chr_change$V1) {
  links_busco_prel[links_busco_prel$seq_id==i, "seq_id"] <- chr_change[chr_change$V1==i, "V2"]
}

unique(links_busco_prel$seq_id)

#only chr number
links_busco_prel$seq_id <- sub(".*_", "", links_busco_prel$seq_id)
#add taxaname to chr number
links_busco_prel$seq_id <- paste(links_busco_prel$bin_id, links_busco_prel$seq_id, sep = "_")

sequences <- unique(links_busco_prel[,c("bin_id", "seq_id", "length")])

#save renamed
write.table(sequences, paste("../intermediate/", prefix, "_sequences_renamed.table", sep = ""))

#save table, for rearrangement stats and syntny analysis
write.table(links_busco_prel, file = paste("../intermediate/", prefix, "_busco_full_table_formated.table", sep = ""))

