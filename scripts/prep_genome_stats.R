#To make a table of current genome stats

#Add variables
#Genome info from concatenated fai (from samtools faidx) 
input_chr_length <- "../input/chr_length_all.tsv"


#read file
chr_length <- read.csv(input_chr_length, sep = "\t", header = F)
colnames(chr_length) <- c("seq_id", "length", "bin_id")

#make taxa (bin) and adjust names to match busco
chr_length$bin_id <-  sub(".*/", "", chr_length$bin_id)
chr_length$bin_id <-  sub("\\.fai", "", chr_length$bin_id)
chr_length$bin_id <-  sub("\\.gz", "", chr_length$bin_id)
unique(chr_length$bin_id)


assembly_stat.df <- data.frame(taxa=0, total_size=0, in_chr=0)

for (i in unique(chr_length$bin_id)) {

assembly_stat.df <- rbind(assembly_stat.df,
                          data.frame(taxa=c(i), 
                               chr_length %>% filter(bin_id==i) %>%
                                 summarise(total_size=sum(length)),
                               chr_length %>% filter(bin_id==i) %>%
                                 filter(str_detect(pattern = "SUPER", seq_id)) %>%
                                 summarise(in_chr=sum(length))))
}

assembly_stat.df$percent <- assembly_stat.df$in_chr/assembly_stat.df$total_size

for (i in unique(chr_length$bin_id)) {
 print( 
chr_length %>% filter(bin_id==i) %>%
  filter(str_detect(pattern = "SUPER", seq_id)) 
)
}
