
# Genomic differentiation between Mech. nesaea compared to Mech. polymnia and Mech. lysimnia

# Read in fd, FST/Dxy and PBS values
fdB<-read.csv("D:/Dropbox/Ithomiines/Mechanitis/nesaea/Mechanitis.nesaea.filtered.2000.fd.polBr_nes_lysBr_mess.csv",sep=",",header=T,na.strings="nan")
fst<-read.csv("D:/Dropbox/Ithomiines/Mechanitis/nesaea/Mechanitis.nesaea.filtered.2000.Fst.Dxy.pi.csv",na.strings = "nan")
pbs<-read.table("D:/Dropbox/Ithomiines/Mechanitis/nesaea/polW_polB_nes.pbs",header=T)

#Chromosomal breakpoints
br<-read.csv("D:/Dropbox/Ithomiines/Mechanitis/nesaea/breakpoints_polymnia_vs_all.csv")

# Remove lines without results 
fdB<-fdB[!is.na(fdB$fdM),]
fst<-fst[!is.na(fst$Fst_M_nesaea_Brazil_M_polymnia_Brazil),]
pbs<-pbs[!is.na(pbs$Chromo),]

# Replace chromosome names to just show the chromosome numbers
pbs$Chromo<-as.integer(sub(pbs$Chromo,pattern = "SUPER_",replacement = ""))
fdB$scaffold<-as.integer(sub(fdB$scaffold,pattern = "SUPER_",replacement = ""))
fst$scaffold<-as.integer(sub(fst$scaffold,pattern = "SUPER_",replacement = ""))

# Add additive columns to plot the chromosomes next to each other
chr<-read.table("D:/Dropbox/Ithomiines/Mechanitis/chromosomes.list",header=T)
rownames(chr)<-chr$chr
chr$mid<-cumsum(chr$length)-chr$length/2

fdB$add<-fdB$mid+chr$add[fdB$scaffold]
fst$add<-fst$mid+chr$add[fst$scaffold]
pbs$add<-pbs$Middle+chr$add[pbs$Chromo]
br$add<-br$br_start+chr$add[br$chr]

# Combining Fst and fd
fst<-merge(fst,fdB[,c("scaffold","start","fdM")],by = c("scaffold","start"))

# introgression outliers:
intr<-fdB[fdB$fdM>quantile(fdB$fdM,0.995),]
intr$start<-intr$start+chr$add[intr$scaffold]
intr$end<-intr$end+chr$add[intr$scaffold]

# Adding information if PBS window (midpoint) is within an introgression window
pbs$intr<-lapply(pbs$add,FUN = function(x) {length(which(with(intr, start<x & intr$end>x)))>0})


##### Figure 3  #####
png("D:/Dropbox/Ithomiines/Mechanitis/nesaea/Fig3_GenomeScans.png",
    width = 9590,height=3290,res=1000)

par(mfrow=c(3,1))
par(xaxs="i",mar=c(0,4,0,1))

# FST nesaea polymnia
plot(fst$add,fst$Fst_M_nesaea_Brazil_M_polymnia_Brazil,col=rep(c("black","darkgrey"),13)[fst$scaffold],
     pch=19,cex=0.5,xlab=NA,ylab=expression('F'['ST']*' nes - pol'),xlim=c(0,257355748),xaxt='n',ylim=c(0,0.9),yaxt="n")
abline(v=intr$add,col="orange",lwd=1.5)
axis(2,at=c(0,0.2,0.4,0.6,0.8),labels = c(0,0.2,0.4,0.6,0.8),las=2)
points(fst$add,fst$Fst_M_nesaea_Brazil_M_polymnia_Brazil,col=rep(c("black","darkgrey"),13)[fst$scaffold],
       pch=19,cex=0.5)
points(fst$add,fst$Fst_M_nesaea_Brazil_M_polymnia_Brazil,
       col=ifelse(fst$Fst_M_nesaea_Brazil_M_polymnia_Brazil>quantile(fst$Fst_M_nesaea_Brazil_M_polymnia_Brazil,0.99),"red",NA),pch=19,cex=1)

# Add locations of chromosomal breakpoints
abline(h=0.8)
rect(xleft = br$add,xright = br$add+(br$br_end-br$br_start),border = NA,
     ybottom=0.8,ytop=1,col="blue")


# PBS nesaea
plot(pbs$add,pbs$PBS_nes_fst,col=rep(c("black","darkgrey"),13)[pbs$Chromo],xlim=c(0,257355748),
     pch=19,cex=0.5,xlab=NA,ylab="PBS Mec. nesaea",xaxt='n',yaxt="n")
abline(v=intr$add,col="orange",lwd=1.5)
axis(2,at=c(0,0.5,1,1.5,2,2.5),labels = c(0,0.5,1,1.5,2,2.5),las=2)
points(pbs$add,pbs$PBS_nes_fst,col=rep(c("black","darkgrey"),13)[pbs$Chromo],
       pch=19,cex=0.5)
points(pbs$add,pbs$PBS_nes_fst,col=ifelse(pbs$PBS_nes_fst>quantile(pbs$PBS_nes_fst,0.99),"red",NA),pch=19,cex=1)

# fdM polymnia Brazil, nesaea, lysimnia Brazil, messenoides
plot(fdB$add,fdB$fdM,col=rep(c("black","darkgrey"),13)[fdB$scaffold],
     pch=19,cex=0.5,xlab=NA,xaxt='n',ylim=c(-0.4,0.8),xlim=c(0,257355748),
     ylab=expression('introgression (f'['dM']*')'),yaxt="n")
abline(v=intr$add,col="orange",lwd=1.5)
axis(2,at=c(-0.2,0,0.2,0.4,0.6,0.8),labels = c(-0.2,0,0.2,0.4,0.6,0.8),las=2)
points(fdB$add,fdB$fdM,col=rep(c("black","darkgrey"),13)[fdB$scaffold],pch=19,cex=0.5,xlab=NA,ylab="fdM_Brazil",xaxt='n',ylim=c(-0.3,1))
points(fdB$add,fdB$fdM,col=ifelse(fdB$fdM>quantile(fdB$fdM,0.99),"red",NA),cex=1,pch=19)

dev.off()


#### Supplementary Figure ####

# Compare windows that are introgression outliers versus the rest for the supplementary figure
png("D:/Dropbox/Ithomiines/Mechanitis/nesaea/SupplementaryFigure_vioplots.png",
    width = 959,height=329,res=100)
par(mfrow=c(1,4),mar=c(4,4,1,1))
require("vioplot")

vioplot(fst$Fst_M_nesaea_Brazil_M_polymnia_Brazil[fst$fdM<quantile(fdB$fdM,0.99)],
        fst$Fst_M_nesaea_Brazil_M_polymnia_Brazil[fst$fdM>quantile(fdB$fdM,0.99)],
        ylab=expression('F'['ST']*' '*italic('Mec. nesaea')*' vs '*italic('Mec. polymnia')),col=c("grey","orange"),
        names=c("normal",expression('top 1% f'['dM'])))
abline(h=median(fst$Fst_M_nesaea_Brazil_M_polymnia_Brazil[fst$fdM>quantile(fdB$fdM,0.99)]))
kruskal.test(fst$Fst_M_nesaea_Brazil_M_polymnia_Brazil,g = fst$fdM>quantile(fdB$fdM,0.99))
mtext("A",las=2,side = 2,line = 3,at=max(fst$Fst_M_nesaea_Brazil_M_polymnia_Brazil))

vioplot(pbs$PBS_nes_fst[pbs$intr==FALSE],
        pbs$PBS_nes_fst[pbs$intr==TRUE],names=c("normal",expression('top 1% f'['dM'])),
        ylab=expression('PBS F'['ST']*' '*italic('Mec. nesaea')),col=c("grey","orange"))
abline(h=median(pbs$PBS_nes_fst[pbs$intr==TRUE]))
mtext("B",las=2,side = 2,line = 3,at=max(pbs$PBS_nes_fst))

vioplot(pbs$PBS_nes_dxy[pbs$intr==FALSE],
        pbs$PBS_nes_dxy[pbs$intr==TRUE],names=c("normal",expression('top 1% f'['dM'])),
        ylab=expression('PBS D'['XY']*' '*italic('Mec. nesaea')),col=c("grey","orange"))
abline(h=median(pbs$PBS_nes_dxy[pbs$intr==TRUE]))
mtext("C",las=2,side = 2,line = 3,at=max(pbs$PBS_nes_dxy))

# nesaea - lysimnia vs nesaea - polymnia
plot(fst$dxy_M_lysimnia_Brazil_M_nesaea_Brazil,fst$dxy_M_polymnia_West_M_nesaea_Brazil,
     col=ifelse(fst$fdM>0.4,"orange","grey"),pch=19,
     xlab=expression('D'['XY']*' '*italic('Mec. nesaea')*' vs '*italic('Mec. lysimnia')),
     ylab=expression('D'['XY']*' '*italic('Mec. nesaea')*' vs '*italic('Mec. polymnia')))
points(fst$dxy_M_lysimnia_Brazil_M_nesaea_Brazil,fst$dxy_M_polymnia_West_M_nesaea_Brazil,
       col=ifelse(fst$fdM>0.4,"darkorange",NA),pch=19,cex=1.5)
mtext("D",las=2,side = 2,line = 3,at=max(fst$dxy_M_polymnia_West_M_nesaea_Brazil,na.rm=T))

dev.off()
