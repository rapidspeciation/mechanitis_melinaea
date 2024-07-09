#### SCRIPT 1 = filtering and converting to phylip, calculating tree-files per window (array-script per chromosome) #####

#!/bin/bash

# (Note: I did this per chromosome, so this was run with an array-script where NR goes through a list of the chromosomes)

NR=`sed -n -e "$LSB_JOBINDEX p" chrom.list` # where $LSB_JOBINDEX comes from submitting to the cluster - this line will be different on another cluster

# input-files 

VCF_IN=/path/vcf/base/Mechanitis.rmINDEL.max0.5N.minDP3.minmeanDP5.maxmeanDP30.rmINDlowD.$NR.vcf.gz
VCF_OUT=/path/filtered/Window.Mech.base.minGQ10.PhyloInd.$NR.vcf.gz
KEEP=/path/scripts/Filtering/IDlists_for_filtering/IDforPhylo
OUT=${VCF_OUT%%.vcf.gz}
ALN_FILE=${OUT}.phylip
TREE_FILE=/path/Phylo.Mech.base.thin500.minGQ10.PhyloInd.ALL.GTR.wB.treefile

PARTITION_FILE=/path/scripts/WindowTree/partition/mechanitis-$NR
## partition file: a RAxML style file per chromosome with all the windows/loci to run for the analysis 
## that looks like this: 
# DNA, part1 = 1-20001
# DNA, part2 = 200001-220001


# Filter with vcftools
module load vcftools/0.1.16--pl5321hd03093a_7

vcftools --gzvcf ${VCF_IN} --minGQ 10 --keep ${KEEP} --recode --stdout | gzip -c > ${VCF_OUT}

# Turn into phylip-file
vcf2phylip.py -i $VCF_OUT -o $OUT.phylip

# infer the locus trees
iqtree2 -s $ALN_FILE -S $PARTITION_FILE -T 4 --prefix /path/Results/WindowTree/$NR.loci

# compute gene concordance factors (per chromosome)
iqtree2 -t $TREE_FILE --gcf /path/Results/WindowTree/$NR.loci.treefile --prefix /path/Results/WindowTree/$NR.concord


##### SCRIPT 2 = obtaining concordance factor in whole genome and TWISST results (no longer array-script per chromosome) #####

#!/bin/bash

### NOTE #### 
## this only works after SCRIPT 1 has already run - it needs a concatenation of all the treefiles as obtained in the 'infer the locus trees'-step
# concatenation with something like: 
# for i in /path/Results/WindowTree/*.loci.treefile ; do cat $i >> /path/Results/WindowTree/all.loci2.treefile ; done
# and then rename the all.loci2.treefile to all.loci.treefile

# input-files
TREE_FILE=/path/Phylo.Mech.base.thin500.minGQ10.PhyloInd.ALL.GTR.wB.treefile
ALL_TREEFILES=/path/Results/WindowTree/all.loci.treefile

# compute gene concordance factors (whole genome)
iqtree2 -t $TREE_FILE --gcf $ALL_TREEFILES --prefix /path/Results/WindowTree/all.loci.concord

# TWISST-analysis

python twisst.py -t $ALL_TREEFILES -w Mechanitis.twisst.lysB-lysE-poly-out.output.weights.csv.gz -g lysimniaBraz -g lysimniaEcu -g polymnia -g Outgroup --method complete \
--groupsFile /path/scripts/TWISST/groups.all.pol1species.short.tsv

python twisst.py -t $ALL_TREEFILES -w Mechanitis.twisst.nes-polyB-lysB-out.output.weights.csv.gz -g nesaea -g polymniaBraz -g lysimniaBraz -g Outgroup --method complete \ 
--groupsFile /path/scripts/TWISST/groups.all.short.tsv

# ETCETERA for all the TWISST-comparisons of interest
# groupsfile is a tab-separated file that has all the individuals and which group/population they are a part of. 
# You don't need to include all the groups in the groupsfile: -g indicates which groups you want to run the analysis for
