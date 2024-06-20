### Introgression statistics

# The first part was run as an array script, running each chromosome in a separate job in parallel.
chr=`sed -n -e "$LSB_JOBINDEX  p" ../scripts/chr.list`

file=Mechanitis.nesaea.filtered.$chr

module load ISG/python/3.7.4

# compute introgression statistics (fd, fdM, ABBA, BABA, etc) with the script by Simon Martin
ABBABABAwindows.py \
-g ../geno/$file.geno.gz -f phased -w 20000 \
 -o $file.2000.fd.polBr_nes_lysBr_mess.csv.gz \
 -m 100 -s 10000 -P1 polymnia_Brazil -P2 nesaea \
 -P3 lysimnia_Brazil -O messenoides \
 -T 4 --minData 0.5 --popsFile nesaea.popfile --writeFailedWindows


# After running all chromosomes, combine the files:
prefix=Mechanitis.nesaea.filtered
suffix=2000.fd.polBr_nes_lysBr_mess

cat <(zcat $prefix.CHR1.$suffix.csv.gz | head -1) \
 <(zcat $prefix.CHR*.$suffix.csv.gz | grep -v scaff) > $prefix.$suffix.csv 


### PBS (Population Branch Statistic)

module load vcftools/0.1.16--pl5321hd03093a_7

chr=`sed -n -e "$LSB_JOBINDEX  p" ../scripts/chr.list`

# extract only the relevant individuals (polymnia and nesaea)
vcftools --vcf ../vcf/Mechanitis.nesaea.filtered.SNPsonly.$chr.renamed.noOut.vcf  \
--keep polymnia_West.txt --keep polymnia_Brazil.txt --keep nesaea.txt \
--recode --out tmp.$chr

# Run PBScan (note, output names must be very short or they are cut off)
pbscan -vcfp tmp.$chr.recode.vcf \
 -pop1 polymnia_West.txt -pop2 polymnia_Brazil.txt -pop3 nesaea.txt  \
 -out polW_polB_nes_${chr} -win 50 -step 51 -div 1 -min 5 -maf 0.05

# Combine the chromosomes
cat <(head -1 polW_polB_nes_CHR1.pbs) <(cat polW_polB_nes_CHR{1..13}.pbs | grep -v Chromo) <(grep -v Chromo polW_polB_nes_Z.pbs) <(grep -v Chromo polW_polB_nes_W.pbs) > dxy_polW_polB_nes.pbs

# pbscan was also run with div 2 to get Dxy values and the PBS FST and DXY versions were then merged:
paste dxy_polW_polB_nes.pbs <(cut -f 6-8 fst_polW_polB_nes.pbs) > polW_polB_nes.pbs

