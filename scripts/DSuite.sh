#!/bin/bash

# input-files
VCF_IN=/path/vcf/base/Mechanitis.rmINDEL.max0.5N.minDP3.minmeanDP5.maxmeanDP30.rmINDlowD.ALL.vcf.gz
VCF_OUT=/path/filtered/Dsuite.Mech.base.minGQ10.minQ10.DSuiteInd.mac2.var.thin100.ALL.vcf.gz
KEEP=/path/scripts/Filtering/IDlists_for_filtering/IDforDsuite

# Filter with vcftools
module load vcftools/0.1.16--pl5321hd03093a_7

vcftools --gzvcf ${VCF_IN} --min-alleles 2 --max-alleles 2 --mac 2 --minGQ 10 --minQ 10 --thin 100 --keep ${KEEP} --recode --stdout | gzip -c > ${VCF_OUT}

# Run Dtrios (Dsuite_tree.nwk is a newick species-tree based on our nuclear phylogenetic tree; 
#		List_Dsuite.txt is a tab-separated text-file with all the individuals in the vcf-file, with the population they are a part of)

Dsuite Dtrios -t /path/scripts/Dsuite/Dsuite_tree.nwk ${VCF_OUT} /path/scripts/Dsuite/List_Dsuite.txt

# Run fbranch
Dsuite Fbranch /path/scripts/Dsuite/Dsuite_tree.nwk /path/scripts/Dsuite/List_Dsuite_tree.txt -Z > /path/Results/Dsuite/fbranch.Z.txt

module load ISG/R/4.1.0

# Set fbranch values to zero if they are not significant
Rscript /path/scripts/Dsuite/removeNonsignDsuite.r -z 3 -i /path/Results/Dsuite/fbranch.Z.txt -o /path/Results/Dsuite/fbranch.signOnly.txt
