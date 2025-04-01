## Load libraries
require(ggplot2)
require(tidyverse)
require(cowplot)
require(tidyverse)
require(ggthemes)



## General settings ************************************************************
maxTrees = 25

cols = c(
  "#0075DC","#e31a1c","#2BCE48","#0075DC","#0075DC",
  "#e31a1c","#0075DC","#e31a1c","#2BCE48","#2BCE48",
  "#e31a1c","#2BCE48","#0075DC","#0075DC","#0075DC",
  "#2BCE48","#0075DC","#0075DC","#e31a1c","#0075DC",
  "#2BCE48","#e31a1c","#2BCE48","#2BCE48","#e31a1c"
)  

# cols = c(
#   "#0075DC","#94FFB5","#FFA405",
#   "#F9E7AC","#FF5005","#005C31",
#   "#00998F","#9DCC00","#2BCE48",
#   "#F3F3F3","#fdbf6f","#b15928",
#   "#FFB7B2","#C20088","#D1A593",
#   "#00998F","#9DCC00","#2BCE48",
#   "#003380","#F0A3FF","#740AFF",
#   "#426600","#C20088","#FF0010",
#   "#FFB7B2","#D1A593","#F3F3F3"
# )  

# cols = c(
#   "#0075DC","#FFA405","#b2df8a",
#   "#33a02c","#fb9a99","#e31a1c",
#   "grey80","grey80","grey80",
#   "grey80","grey80","grey80",
#   "grey80","grey80","grey80"
# )  


## Read BPP results ************************************************************
dir = "/n/holyscratch01/mallet_lab/fseixas/ithomiini/mechanitis/mec_allspp/bpp-a01-noncod.l300_g2000/"
fil = "bpp_block_treeCounts.txt"
bpp = read.table(paste0(dir,fil))
colnames(bpp) = c("chromRegion","blockNb","replicate","tree","count")
# Add info
bpp$chrom = str_replace(bpp$chromRegion, "SUPER_", "")
bpp$Region = str_replace(bpp$chromRegion, "SUPER_", "")
bpp$name = paste0(bpp$chromRegion,".",bpp$blockNb,".",bpp$replicate) 
bpp$blockCode = paste0(bpp$chromRegion,".",bpp$blockNb)


## Define block positions in chromosomes ****************************************
tmpbpp = data.frame()
for (ch in unique(bpp$chrom)) {
  s1 = subset(bpp, chrom == ch)
  s1$blockPos = NA
  unq_blocks = unique(s1$blockCode)
  df = data.frame(
    blockCode = unq_blocks, 
    blockPos = seq(1,length(unq_blocks),1)
  )
  for (rw in 1:nrow(df)) {
    bl = df$blockCode[rw]
    nb = df$blockPos[rw]
    s1$blockPos = ifelse(s1$blockCode == bl, nb, s1$blockPos)
  }
  tmpbpp = rbind(tmpbpp, s1)
}
bpp = tmpbpp

## Get unique trees and count the Proportion of each three *********************
uniqueTrees = unique(bpp$tree)
GlobalCount = data.frame(Tree = uniqueTrees, Count = 0)
for (rw in 1:nrow(GlobalCount)) {
  tr = GlobalCount$Tree[rw]
  s1 = subset(bpp, tree == tr)
  ct = sum(s1$count)
  GlobalCount$Count[rw] = ct
}
GlobalCount$Proportion = GlobalCount$Count / sum(GlobalCount$Count)
GlobalCount = GlobalCount[order(GlobalCount$Proportion, decreasing = T),]
maxTrees = min(nrow(GlobalCount), maxTrees)
GlobalCountMajor = GlobalCount[1:maxTrees,]
MajorTrees = GlobalCountMajor$Tree
MajorTrees = unique(MajorTrees)
cols = cols[1:length(MajorTrees)]

# write major topologies to file
write.table(x = paste0(MajorTrees,";"), file = paste0(dir,"major_topos.txt"), row.names=F, col.names = F, quote=F)

## MAP trees *******************************************************************
MAPtrees = data.frame()
for (blockRep in unique(bpp$blockCode)) {
  blockRepTrees = subset(bpp, blockCode == blockRep)
  blockRepTrees$Proportion = blockRepTrees$count/sum(blockRepTrees$count)
  MAP = as.character(subset(blockRepTrees, Proportion == max(blockRepTrees$Proportion))$tree)
  MAPprop = subset(blockRepTrees, Proportion == max(blockRepTrees$Proportion))$Proportion*100
  map_df = data.frame(block = blockRep, tree = MAP, Proportion = MAPprop)
  MAPtrees = rbind(MAPtrees, map_df)
}
# Summary of MAP trees
MAPtrees_Summary = as.data.frame(table(MAPtrees$tree))
MAPtrees_Summary = MAPtrees_Summary[order(MAPtrees_Summary$Freq, decreasing = T),]
colnames(MAPtrees_Summary) = c("tree","Frequency")

# recode trees depending on being major or minor
bpp$treeRecode = bpp$tree
bpp$treeRecode = ifelse(bpp$tree %in% MajorTrees, bpp$treeRecode, "Other")
# order trees
bpp$tree = factor(bpp$tree, levels = GlobalCount$Tree)
bpp$treeRecode = factor(bpp$treeRecode, levels = c(MajorTrees,"Other"))
bpp$Proportion = bpp$count / 100000



## Plot Major topologies *******************************************************
library(ape)
library(ggtree)

# read topologies file
topos = read.tree(file=paste0(dir,"major_topos.txt"))

# plot major trees
for (n in 1:length(topos)) {
  ## plot
  tp = write.tree(topos[[n]], digits = 10) # print topology as newick
  tp = gsub(";", "", tp)
  pp = subset(GlobalCount, Tree == tp[1] )
  prop = paste0(sprintf("%.2f", pp$Proportion*100),"%")
  p = ggtree(topos[n], layout = "slanted", col=cols[n], lwd=1) +
    geom_tiplab(size=4) +
    xlim(0,12) +
    ggtitle(paste0("topo",n,": (",prop,")")) +
    theme(text = element_text(size=8))
  assign(paste0("pl",n), p)
}


## Plot BPP along chromosomes **************************************************
bpp$chrom = factor(bpp$chrom, c(seq(1,14,1),"Z"))
# bpp$tree = factor(bpp$tree, levels = MajorTrees)
bpp$treeRecode = factor(bpp$treeRecode, levels = c(MajorTrees,"Other"))
l_plot = bpp %>%
  # filter(chrom == 9) %>% 
  ggplot(aes(fill=treeRecode, y=count, x=blockPos)) + 
  geom_col(position = "stack", col=alpha("black",0.0), width = 1) +
  # geom_bar(aes(group=1), position="stack", stat="identity", alpha=.8, col=alpha("black",0.0), width = 0.95) +
  facet_wrap(~chrom, strip.position="left", ncol=1) +
  scale_x_continuous(expand=c(0,0), breaks=seq(0,300,10), name="Block Number") +
  scale_y_continuous(name="Count") +
  scale_fill_manual(values = c(cols,"grey") ) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text.y = element_blank(),
    axis.line = element_blank(), 
    axis.ticks = element_blank(),
    strip.text.y.left = element_text(angle = 0), 
    strip.background = element_blank(),
    strip.placement = "inside"
  )
r_plot = plot_grid(
  pl1,pl2,pl3,pl4,pl5,
  pl6,pl7,pl8,pl9,pl10,
  pl11,pl12,pl13,pl14,pl15,
  pl16,pl17,pl18,pl19,pl20,
  pl21,pl22,pl23,pl24,pl25,
  nrow=5
)
plot_grid(l_plot,r_plot, rel_widths = c(3,3))






## Plot by chromosome **********************************************************
l_plot = bpp %>%
  group_by(chrom, treeRecode) %>%
  summarise(meanProp = mean(Proportion)) %>%
  # plot
  ggplot(aes(fill=treeRecode, y=meanProp, x=chrom)) + 
  geom_bar(
    aes(group=1), 
    position="fill", stat="identity", alpha=.8, col="white", width = 0.95, lwd=.5
  ) +
  scale_x_discrete(expand=c(0,0), name="Chromosome") +
  scale_y_continuous(name="Count", expand=c(0,0)) +
  # scale_fill_tableau(palette = "Tableau 20") +
  scale_fill_manual(values = c(cols,"grey") ) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle=90, hjust=0, vjust=0.5),
    axis.text.y = element_blank(),
    axis.line = element_blank(), 
    axis.ticks = element_blank(),
    strip.text.y.left = element_text(angle = 0), 
    strip.background = element_blank(),
    strip.placement = "inside"
  )
r_plot = plot_grid(
  pl1,pl2,pl3,pl4,pl5,
  pl6,pl7,pl8,pl9,pl10,
  pl11,pl12,pl13,pl14,pl15,
  nrow=3
)
plot_grid(l_plot,r_plot, rel_widths = c(2,3))




