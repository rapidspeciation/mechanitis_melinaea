## Load libraries
require(ggplot2)
require(tidyverse)
require(cowplot)
require(tidyverse)
require(ggthemes)
library(ape)
library(ggtree)

dir = d
maxTrees = 25
cols = allCols


## Functions
plotBPP = function(dir, maxTrees, cols, topLevels) {
  
  ## Read BPP results ************************************************************
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


  # name trees
  bpp$treeName = "Other"
  for (n in 1:length(MajorTrees)) {
    bpp$treeName = ifelse(bpp$treeRecode == MajorTrees[n], paste0("topo",n), bpp$treeName)
  }
  bpp$treeName = factor(bpp$treeName, levels = topLevels)
  
  
  ## Plot Major topologies *******************************************************
  
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
      theme(text = element_text(size=9))
    assign(paste0("pl",n), p)
  }
  
  ## Plot BPP along chromosomes **************************************************
  # bpp$chrom = factor(bpp$chrom, seq(1,14,1))
  bpp$chrom = factor(bpp$chrom, c(seq(1,14,1),"Z"))
  bpp$tree = factor(bpp$tree, levels = MajorTrees)
  bpp$treeRecode = factor(bpp$treeRecode, levels = c(MajorTrees,"Other"))
  
  max_x = floor(max(bpp$blockNb)/10)*10
  
  
  l_plot = bpp %>%
    # ggplot(aes(fill=treeRecode, y=count, x=blockPos)) +
    ggplot(aes(fill=treeName, y=count, x=blockPos)) +
    geom_vline(xintercept = seq(0,max_x,10), lty=3) +
    geom_bar(
      position="fill", stat="identity", alpha=.8, , col=alpha("black",0.0), width = 1.0, lwd=.0
    ) +
    # geom_col(position = "stack", col=alpha("black",0.0), width = 1) +
    facet_wrap(~chrom, strip.position="left", ncol=1) +
    scale_x_continuous(expand=c(0,0), breaks=seq(0,300,10), name="Block Number") +
    scale_y_continuous(name="Count") +
    scale_fill_manual(values = c(cols,"grey90") ) +
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

  r1_plot = plot_grid(
    pl1,pl2,pl3,pl4,pl5,
    pl6,pl7,pl8,pl9,pl10,
    pl11,pl12,pl13,pl14,pl15,
    pl16,pl17,pl18,pl19,pl20,
    pl21,pl22,pl23,pl24,pl25,
    nrow=5
  )
  
  r2_plot = bpp %>%
    group_by(chrom, treeName) %>%
    # summarise(meanProp = mean(Proportion)) %>%
    summarise(sumProp = sum(Proportion)) %>%
    # plot
    ggplot(aes(fill=treeName, y=sumProp, x=chrom)) + 
    geom_bar(
      aes(group=1), 
      position="fill", stat="identity", alpha=.8, col="white", width = 0.95, lwd=.0
    ) +
    scale_x_discrete(expand=c(0,0), name="Chromosome") +
    scale_y_continuous(name="Count", expand=c(0,0)) +
    scale_fill_manual(values = c(cols,"grey90") ) +
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
  
  rig = plot_grid(r1_plot,r2_plot, nrow = 2, rel_heights = c(3,2))
  p = plot_grid(l_plot, rig, ncol=2)
  return(p)
}


## Define colors
mel_outgroup = c(
  "#0075DC", "#FFA405", "#FFA405", "#0075DC", "#0075DC",
  "#0075DC", "#FFA405", "#FFA405", "#0075DC", "#0075DC",
  "#FFA405", "#b2df8a", "#0075DC", "#FFA405", "#FFA405",
  "#FFA405", "#b2df8a", "#fb9a99", "#D1A593", "#0075DC",
  "#FFA405", "#b2df8a", "#0075DC", "#0075DC", "#FFA405"
)

mel_outgroup2 = c(
  "#0075DC", "#FFA405", "#0075DC", "#FFA405", "#0075DC",
  "#FFA405", "#0075DC", "#FFB7B2", "#0075DC", "#FFA405",
  "#0075DC", "#FFA405", "#0075DC", "#FFA405", "#b2df8a",
  "#b2df8a", "#FFA405", "#0075DC", "#FFA405", "#FFA405",
  "#b2df8a", "#b2df8a", "#0075DC", "#b2df8a", "#0075DC"
)



mec_allspp = c(
   "topo1"="#006d2c", "topo2"="#006d2c", "topo3"="#006d2c", "topo4"="#993404",  "topo5"="#006d2c",
   "topo6"="#993404", "topo7"="#1d91c0", "topo8"="#41ae76", "topo9"="#993404", "topo10"="#41ae76",
  "topo11"="#1d91c0","topo12"="#1d91c0","topo13"="#fdae6b","topo14"="#41ae76", "topo15"="#99d8c9",
  "topo16"="#fdae6b","topo17"="#253494","topo18"="#f16913","topo19"="#fdae6b", "topo20"="#FFD700",
  "topo21"="#FFD700","topo22"="#253494","topo23"="#99d8c9","topo24"="#253494", "topo25"="#99d8c9", 
  "Other"="grey90"
)

## mechanitis colors
mec_allspp = c(
  "#006d2c","#006d2c","#006d2c","#993404","#006d2c",
  "#993404","#1d91c0","#41ae76","#993404","#41ae76",
  "#1d91c0","#1d91c0","#fdae6b","#41ae76","#99d8c9",
  "#fdae6b","#253494","#f16913","#fdae6b","#FFD700",
  "#FFD700","#253494","#99d8c9","#253494","#99d8c9"
)

mec_allspp = c(
   "topo1"="#006d2c", "topo2"="#006d2c", "topo3"="#006d2c", "topo4"="#993404",  "topo5"="#006d2c",
   "topo6"="#993404", "topo7"="#1d91c0", "topo8"="#41ae76", "topo9"="#993404", "topo10"="#41ae76",
  "topo11"="#1d91c0","topo12"="#1d91c0","topo13"="#fdae6b","topo14"="#41ae76", "topo15"="#99d8c9",
  "topo16"="#fdae6b","topo17"="#253494","topo18"="#f16913","topo19"="#fdae6b", "topo20"="#FFD700",
  "topo21"="#FFD700","topo22"="#253494","topo23"="#99d8c9","topo24"="#253494", "topo25"="#99d8c9", 
  "Other"="grey90"
)


mec_allspp2 = c(
  "#0075DC","#FFA405","#b2df8a","#0075DC","#0075DC",
  "#FFA405","#0075DC","#FFA405","#b2df8a","#b2df8a",
  "#FFA405","#b2df8a","#0075DC","#0075DC","#0075DC",
  "#b2df8a","#0075DC","#0075DC","#FFA405","#0075DC",
  "#b2df8a","#FFA405","#b2df8a","#b2df8a","#FFA405"
)

lis_nes_pol = c(
  "#0075DC","#FFA405","#e31a1c","#0075DC","#fb9a99",
  "#0075DC","#0075DC","#ff7f00","#0075DC","#6a3d9a",
  "#F9E7AC","#505050","#D62728","#9467BD","#D1A593",
  "#00998F","#9DCC00","#2BCE48","#003380","#F0A3FF",
  "#740AFF","#426600","#C20088","#FF0010","#FFB7B2"
)

allCols = c(
  "#1F77B4", "#FF7F0E", "#2CA02C", "#D62728", "#9467BD",
  "#8C564B", "#E377C2", "#7F7F7F", "#BCBD22", "#17BECF",
  "#AEC7E8", "#FFBB78", "#98DF8A", "#FF9896", "#C5B0D5",
  "#C49C94", "#F7B6D2", "#C7C7C7", "#DBDB8D", "#9EDAE5",
  "#393B79", "#637939", "#8C6D31", "#843C39", "#7B4173"
)


## Plot BPP  ************************************************************

## lysimnia-nesaea-polymnia
d = "/n/holyscratch01/mallet_lab/fseixas/ithomiini/mechanitis/lis_nes_pol/bpp-a01-noncod.l300_g2000/"
plotBPP(d, 25, allCols)
plotBPP(d, 25, lis_nes_pol)

## mechanatis
d = "/n/holyscratch01/mallet_lab/fseixas/ithomiini/mechanitis/mec_allspp/bpp-a01-noncod.l300_g2000/"
plotBPP(d, 25, allCols)
plotBPP(d, 25, mec_allspp, mecLevels)
plotBPP(d, 25, mec_allspp2)

## melinaea
d = "/n/holyscratch01/mallet_lab/fseixas/ithomiini/melinaea/mel_outgroup/bpp-a01-noncod.l300_g2000/"
plotBPP(d, 25, allCols)
plotBPP(d, 25, mel_outgroup)

d = "/n/holyscratch01/mallet_lab/fseixas/ithomiini/melinaea/mel_outgroup2/bpp-a01-noncod.l300_g2000/"
plotBPP(d, 25, allCols)
plotBPP(d, 25, mel_outgroup2)




