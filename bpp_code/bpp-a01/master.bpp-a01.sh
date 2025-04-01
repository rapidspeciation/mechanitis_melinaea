## Generate Dataset ********************************************************

# Arguments
popset="lis_nes_pol"
lsize=300
lstep=2000
window=`echo "l${lsize}_g${lstep}"`
refnam="mechanitis"
export popset=$popset 
export refnam=$refnam 
export window=$window

mkdir $popset; cd $popset

## 1. prepare bed files
sbatch bpp_code/1.bpp_prepareBed.slurm \
    ${refnam} 2000 \
    250 50 \
    $lsize $lstep \
    /n/holyscratch01/mallet_lab/fseixas/ithomiini/${refnam}/0.data/ \
    ~/code/heliconius_seixas/ithomiini/popset/${popset}.pset

## 2. produce alignments
sbatch bpp_code/2.reg2aln.slurm \
    ${refnam} a01 noncod.${window} \
    0.80 0.0 1.0 30 \
    ~/code/heliconius_seixas/ithomiini/popset/${popset}.imap \
    ~/code/heliconius_seixas/ithomiini/popset/${popset}.pmin \
    /n/holyscratch01/mallet_lab/fseixas/ithomiini/${refnam}/reference/${refnam}.repeats.txt

mkdir h0.80_v0.0_i1.0_l30; mv 2.vcf2bpp/ 3.0.blocks/ 3.1.blocksAlign/ 3.2.bpp/ h0.80_v0.0_i1.0_l30
mkdir h0.50_v0.0_i0.5_l50; mv 2.vcf2bpp/ 3.0.blocks/ 3.1.blocksAlign/ 3.2.bpp/ h0.50_v0.0_i0.5_l50
mkdir h0.75_v0.2_i0.5_l50; mv 2.vcf2bpp/ 3.0.blocks/ 3.1.blocksAlign/ 3.2.bpp/ h0.75_v0.2_i0.5_l50
mkdir h0.75_v0.2_i0.5_l50; mv 2.vcf2bpp/ 3.0.blocks/ 3.1.blocksAlign/ 3.2.bpp/ h0.75_v0.2_i0.5_l50


## define blocks to analyze
bash bpp_code/bpp-a01/2.locusBlocks.sh \
    ${refnam} \
    noncod \
    ${window} \
    100 \
    none \
    /n/holyscratch01/mallet_lab/fseixas/ithomiini/${refnam}/reference/${refnam}.scaffold.sizes
cd bpp-a01-noncod.${window}/
awk '{sum=sum+$4} END {print sum/NR/1000000}' blockSummary.txt


## BPP control files
cd bpp-a01-noncod.${window}/

# mechanatis
bash bpp_code/3.runBPP.sh \
    "7 lysE mac maz men mes nes polW" \
    "1 1 1 1 1 1 1" \
    "((mac,maz),((polW,men),(lysE,(nes,mes))));" \
    "1 1 1 1 1 1 1" \
    bpp_code/popset/${popset}.imap \
    CAM071978 1

# melinaea all
bash bpp_code/3.runBPP.sh \
    "11 ide iso lil lud mae mar men mne sat tar mot" \
    "1 1 1 1 1 1 1 1 1 1 1 " \
    "(((((mae,mot),(lil,lud)),ide),(mne,((mar,tar),sat))),(iso,men));" \
    "1 1 1 1 1 1 1 1 1 1 1 " \
    bpp_code/popset/${popset}.imap \
    CAM036282 1

# melinea ingroup
bash bpp_code/3.runBPP.sh \
    "7 mar mae sat tar mot men mne" \
    "1 1 1 1 1 1 1 " \
    "((mar, mae), ((sat, tar), (mot, (men, mne))));" \
    "1 1 1 1 1 1 1 " \
    ~/code/heliconius_seixas/ithomiini/popset/${popset}.imap \
    MJ094199 1

# melinea outgroup
bash bpp_code/3.runBPP.sh \
    "7 lud lil iso ide mne mot mar" \
    "1 1 1 1 1 1 1 " \
    "((lud, lil), ((iso, ide), (mne, (mot, mar))));" \
    "1 1 1 1 1 1 1 " \
    ~/code/heliconius_seixas/ithomiini/popset/${popset}.imap \
    MJ094199 1

# melinea outgroup
bash bpp_code/3.runBPP.sh \
    "6 lil iso ide mne mot mar" \
    "1 1 1 1 1 1 " \
    "((lil, iso), ((ide, mne), (mot, mar)));" \
    "1 1 1 1 1 1 " \
    ~/code/heliconius_seixas/ithomiini/popset/${popset}.imap \
    MJ094199 1


# melinea outgroup
bash bpp_code/3.runBPP.sh \
    "5 lysBra lysEcu nes polBra polWes" \
    "1 1 1 1 1 " \
    "((lysBra, lysEcu), (nes, (polBra, polWes)));" \
    "1 1 1 1 1 " \
    ~/code/heliconius_seixas/ithomiini/popset/${popset}.imap \
    CAM071978 1


## launch BPP
ls 3.1.blocksAlign/ | sed 's,.seqs.txt,,'g | egrep -wv "SUPER_Z" | xargs -n 1 sh -c '
    sbatch bpp_code/bpp01/launch_bpp-a01.slurm $0 rep1
'

## check progress
wc -l 3.2.bpp/SUPER_8.block_???.rep1/SUPER_*.block_???.rep1.mcmc.list



## Summary of BPP-A01 *****************************************************************

#/ check if files finished or need to be resumed; if done copy final files to finalized folder
mkdir 3.2.bpp/finished/
bash bpp_code/bpp_status.sh 100001
#/ concatenate output files // tree Counts
rm -f bpp_block_treeCounts.txt
ls 3.2.bpp/finished/*out* | xargs -n 1 -P 1 sh -c '
    perl bpp_code/scripts/parseBPP_output.pl $0 \
    >> bpp_block_treeCounts.txt
'
# transform file to more readable format
sed -i 's,3.2.bpp/finished/,,g' bpp_block_treeCounts.txt
sed -i 's,.out.txt,,g' bpp_block_treeCounts.txt
sed -i 's,\.,\t,g' bpp_block_treeCounts.txt
sed -i 's,block\_,,g' bpp_block_treeCounts.txt
