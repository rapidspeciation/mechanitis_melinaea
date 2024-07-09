##### Script for making two phylogenies and comparing them #####

#### Load libraries ####
require("ape")
require(phytools)
library(colorblindcheck)
library(MoMAColors)

#### Melinaea ####
# choose colour palette that's colour blind friendly
pal <- c("#000000","#004949","#009292","#24ff24",
         "#ff6db6","#ffb6db","#490092","#b66dff",
         "#006ddb","#6db6ff","#920000","#924900","#db6d00")
palette_check(pal, plot = TRUE)
palette_dist(pal, cvd="deu")

# Read in the trees to compare 
mel<-read.tree("PhD/R-scripts/Melanitis/trees/Phylo.Mel.base.thin500.minGQ10.PhyloInd.ALL.GTR.wB.treefile.renamed")
mel_mt<-read.tree("PhD/R-scripts/Melanitis/trees/Phylo.Mel.base.minGQ10.PhyloInd.MT.GTR.wB.treefile.renamed")

mel <- root(mel,c("CAM072004_outgroup_Ecuador-East","FS50864856_outgroup_Ecuador-West"))
mel_mt <- root(mel_mt,c("CAM072004_outgroup_Ecuador-East","FS50864856_outgroup_Ecuador-West"))

mel <- ladderize(mel)
mel_mt <- ladderize(mel_mt)

# Get a dataframe with the association (here just twice the same individual)
assoc<-cbind(mel$tip.label, mel$tip.label)

# Get the cophyloplot object
obj<-cophylo(mel,mel_mt,assoc=assoc)

# Get the tiplabels in the new order
nucTips<-obj$trees[[1]]$tip.label
mtDNATips<-obj$trees[[2]]$tip.label

### make figure before renaming ####
# get species colours
speciesCol<-as.data.frame(assoc[,1])
names(speciesCol)<-"ind"
speciesCol$col<-ifelse(grepl(speciesCol$ind,pattern="idae"),"#ff6db6",
                       ifelse(grepl(speciesCol$ind,pattern="mneme"),"#ffb6db",
                              ifelse(grepl(speciesCol$ind,pattern="isocomma"),"#490092",
                                     ifelse(grepl(speciesCol$ind,pattern="lilis"),"#920000",
                                            ifelse(grepl(speciesCol$ind,pattern="ludovica"),"#6db6ff",
                                                   ifelse(grepl(speciesCol$ind,pattern="satevis"),"#24ff24",
                                                          ifelse(grepl(speciesCol$ind,pattern="marsaeus"),"#006ddb",
                                                                 ifelse(grepl(speciesCol$ind,pattern="menophilus"),"#db6d00",
                                                                        ifelse(grepl(speciesCol$ind,pattern="zaneka"),"#db6d00",
                                                                               ifelse(grepl(speciesCol$ind,pattern="outgroup"),"#000000",
                                                                                      ifelse(grepl(speciesCol$ind,pattern="mothone"),"#b66dff","#24ff24")))))))))))



# Plot a cophyloplot
plot.cophylo(obj,link.type="curved",link.lwd=3,fsize=0.3,pts=F,ftype="off",
             link.lty="solid",link.col=speciesCol$col)

# Add bootstrap support symbols to both trees
#-> bootstrap 100 = large dot, bootstrap >90 is small dot, below 90 no dot.

nodelabels.cophylo(pch=15,
                   col=ifelse(as.integer(obj$trees[[1]]$node.label)==100,"#41424C",
                              ifelse(as.integer(obj$trees[[1]]$node.label)>90,"#41424C",NA)),
                   cex=ifelse(as.integer(obj$trees[[1]]$node.label)==100,0.7,
                              ifelse(as.integer(obj$trees[[1]]$node.label)>90,0.4,NA)))
nodelabels.cophylo(pch=15,col=ifelse(as.integer(obj$trees[[2]]$node.label)==100,"#41424C",
                                     ifelse(as.integer(obj$trees[[2]]$node.label)>90,"#41424C",NA)),
                   cex=ifelse(as.integer(obj$trees[[2]]$node.label)==100,0.7,
                              ifelse(as.integer(obj$trees[[2]]$node.label)>90,0.4,NA)),
                   which="right") 

# nuclear tree
tiplabels.cophylo(pch=16,cex=0.8,which = "left",
                  col=ifelse(grepl(nucTips,pattern="ludovica"),"#6db6ff",
                             ifelse(grepl(nucTips,pattern="lilis"),"#920000",
                                    ifelse(grepl(nucTips,pattern="isocomma"),"#490092",
                                           ifelse(grepl(nucTips,pattern="idae"),"#ff6db6",
                                                  ifelse(grepl(nucTips,pattern="tarapotensis"),"#24ff24",
                                                         ifelse(grepl(nucTips,pattern="satevis"),"#24ff24",
                                                                ifelse(grepl(nucTips,pattern="mothone"),"#b66dff",
                                                                       ifelse(grepl(nucTips,pattern="menophilus"),"#db6d00",
                                                                              ifelse(grepl(nucTips,pattern="zaneka"),"#db6d00",
                                                                                     ifelse(grepl(nucTips,pattern="marsaeus"),"#006ddb",
                                                                                            ifelse(grepl(nucTips,pattern="maeonis"),"#24ff24",
                                                                                                   ifelse(grepl(nucTips,pattern="mneme"),"#ffb6db",
                                                                                                          ifelse(grepl(nucTips,pattern="SAN2500"),"#24ff24","#000000"))))))))))))))


# mitochondrial tree
tiplabels.cophylo(pch=16,cex=0.8,which = "right",
                  col=ifelse(grepl(mtDNATips,pattern="ludovica"),"#6db6ff",
                             ifelse(grepl(mtDNATips,pattern="lilis"),"#920000",
                                    ifelse(grepl(mtDNATips,pattern="isocomma"),"#490092",
                                           ifelse(grepl(mtDNATips,pattern="idae"),"#ff6db6",
                                                  ifelse(grepl(mtDNATips,pattern="tarapotensis"),"#24ff24",
                                                         ifelse(grepl(mtDNATips,pattern="satevis"),"#24ff24",
                                                                ifelse(grepl(mtDNATips,pattern="mothone"),"#b66dff",
                                                                       ifelse(grepl(mtDNATips,pattern="menophilus"),"#db6d00",
                                                                              ifelse(grepl(mtDNATips,pattern="zaneka"),"#db6d00",
                                                                                     ifelse(grepl(mtDNATips,pattern="marsaeus"),"#006ddb",
                                                                                            ifelse(grepl(mtDNATips,pattern="maeonis"),"#24ff24",
                                                                                                   ifelse(grepl(mtDNATips,pattern="mneme"),"#ffb6db",
                                                                                                          ifelse(grepl(mtDNATips,pattern="SAN2500"),"#24ff24","#000000"))))))))))))))






# Add legends
legend("bottomright",legend=c("Eutresis/Olyras","Mel. ludovica","Mel. lilis","Mel. isocomma","Mel. idae",
                              "Mel. satevis","Mel. mothone","Mel. menophilus",
                              "Mel. marsaeus","Mel. mneme"),
       fill=                c("#000000","#6db6ff","#920000","#490092","#ff6db6",
                              "#24ff24","#b66dff","#db6d00",
                              "#006ddb","#ffb6db"),
       bty="n",title="Species",cex=1,xpd=NA)


##### Mechanitis ####

#### Read in the data ####
# PhyloInd
mech_mt<-read.tree("PhD/R-scripts/Melanitis/trees/Phylo.Mech.base.minGQ10.MT.PhyloInd.wB.treefile.renamed")
mech<-read.tree("PhD/R-scripts/Melanitis/trees/Phylo.Mech.base.thin500.minGQ10.PhyloInd.ALL.wB.treefile.renamed")

# Root Mechanitis on Forbestra
mech <- midpoint.root(mech)
mech_mt <- midpoint.root(mech_mt)

# Ladderize the trees
mech<-ladderize(mech)
mech_mt<-ladderize(mech_mt)

# Get a dataframe with the association (here just twice the same individual)
assoc<-cbind(mech$tip.label, mech$tip.label)

obj<-cophylo(mech,mech_mt,assoc=assoc,cex=0.1)

# Get the tiplabels in the new order
nucTips<-obj$trees[[1]]$tip.label
mtDNATips<-obj$trees[[2]]$tip.label

#### Make the cophylo-plot ####
## choose colour palette
display.all.moma(7, colorblind_only = TRUE)
Rat <- moma.colors("Rattner", n = 7)

speciesCol<-as.data.frame(assoc[,1])
names(speciesCol)<-"ind"
speciesCol$col<-ifelse(grepl(speciesCol$ind,pattern="forbestra"),"black",
                                            ifelse(grepl(speciesCol$ind,pattern="macrinus"),Rat[5],
                                                   ifelse(grepl(speciesCol$ind,pattern="mazaeus"),Rat[7],
                                                          ifelse(grepl(speciesCol$ind,pattern="menapis"),Rat[6],
                                                                 ifelse(grepl(speciesCol$ind,pattern="messenoides"),Rat[4],
                                                                        ifelse(grepl(speciesCol$ind,pattern="lysimnia"),Rat[5],
                                                                               ifelse(grepl(speciesCol$ind,pattern="nesaea"),Rat[5],
                                                                                      ifelse(grepl(speciesCol$ind,pattern="polymnia"),Rat[1],"white"))))))))

# Plot a cophyloplot
plot.cophylo(obj,link.type="curved",link.lwd=3,fsize=0.3,pts=F,ftype="off",
             link.lty="solid",link.col=speciesCol$col)

nodelabels.cophylo(pch=15,
                   col=ifelse(as.integer(obj$trees[[1]]$node.label)==100,"#41424C",
                              ifelse(as.integer(obj$trees[[1]]$node.label)>90,"#41424C",NA)),
                   cex=ifelse(as.integer(obj$trees[[1]]$node.label)==100,0.7,
                              ifelse(as.integer(obj$trees[[1]]$node.label)>90,0.4,NA)))
nodelabels.cophylo(pch=15,col=ifelse(as.integer(obj$trees[[2]]$node.label)==100,"#41424C",
                                     ifelse(as.integer(obj$trees[[2]]$node.label)>90,"#41424C",NA)),
                   cex=ifelse(as.integer(obj$trees[[2]]$node.label)==100,0.7,
                              ifelse(as.integer(obj$trees[[2]]$node.label)>90,0.4,NA)),
                   which="right") 

# Add tip labels with species colours
# mitochondrial tree
tiplabels.cophylo(pch=16,cex=0.8,which = "left",
                  col=ifelse(grepl(nucTips,pattern="forbestra"),"black",
                                                 ifelse(grepl(nucTips,pattern="macrinus"),Rat[5],
                                                         ifelse(grepl(nucTips,pattern="mazaeus"),Rat[7],
                                                                ifelse(grepl(nucTips,pattern="menapis"),Rat[6],
                                                                       ifelse(grepl(nucTips,pattern="messenoides"),Rat[4],
                                                                              ifelse(grepl(nucTips,pattern="lysimnia"),Rat[5],
                                                                                     ifelse(grepl(nucTips,pattern="nesaea"),Rat[5],
                                                                                            ifelse(grepl(nucTips,pattern="polymnia"),Rat[1],"white")))))))))

tiplabels.cophylo(pch=16,cex=0.8,which = "right",
                  col=ifelse(grepl(mtDNATips,pattern="forbestra"),"black",
                                                  ifelse(grepl(mtDNATips,pattern="macrinus"),Rat[5],
                                                         ifelse(grepl(mtDNATips,pattern="mazaeus"),Rat[7],
                                                                ifelse(grepl(mtDNATips,pattern="menapis"),Rat[6],
                                                                       ifelse(grepl(mtDNATips,pattern="messenoides"),Rat[4],
                                                                              ifelse(grepl(mtDNATips,pattern="lysimnia"),Rat[5],
                                                                                     ifelse(grepl(mtDNATips,pattern="nesaea"),Rat[5],
                                                                                            ifelse(grepl(mtDNATips,pattern="polymnia"),Rat[1],"white")))))))))

legend("bottomright",legend=c("Forbestra","Mech. mazaeus","Mech. menapis","Mech. messenoides","Mech. lysimnia","Mech. polymnia"),
       fill=c("black",Rat[7],Rat[6],Rat[4],Rat[5],Rat[1]), bty="n",title="Species",cex=1)
