# Some files need to ber created before running this script

#/ create keep, pset, pmap and pmin files 
# keep: list of individuals for the analysis (as necessary in plink format)
# pset: list of individuals for the analysis
# pmap: map individuals to groups
# pmin: min. number of individuals retained per population after filters



#/ Generate newick tree for the msci model
~/software/bpp.v4.6.2/src/bpp \
        --msci-create bpp_code/models/msci.modelA.txt | \
        tail -n 1 > bpp_code/models/msci.modelA.nwk

~/software/bpp.v4.6.2/src/bpp \
        --msci-create bpp_code/models/msci.modelFc.txt | \
        tail -n 1 > bpp_code/models/msci.modelFc.rep1.nwk

~/software/bpp.v4.6.2/src/bpp \
        --msci-create bpp_code/models/msci.modelFb.txt | \
        tail -n 1 > bpp_code/models/msci.modelFb.rep1.nwk

~/software/bpp.v4.6.2/src/bpp \
        --msci-create bpp_code/models/msci.modelB1a.txt | \
        tail -n 1 > bpp_code/models/msci.modelB1a.rep1.nwk


## Generate Dataset ********************************************************

# Arguments
popset="ilm_msci"; refnam="melinaea"
# popset="mmp_msci"; refnam="mechanitis"
# popset="1indEcu"; refnam="mechanitis"

lsize=2000
lstep=10000
window=`echo "l${lsize}_g${lstep}"`

export popset=$popset 
export refnam=$refnam 
export window=$window 

mkdir $popset; cd $popset

## 1. prepare bed files
sbatch bpp_code/1.bpp_prepareBed.slurm \
    ${refnam} 5000 \
    250 50 \
    $lsize $lstep \
    ${refnam}/0.data/ \
    bpp_code/popset/${popset}.pset

## 2. produce alignments 
sbatch bpp_code/2.reg2aln.slurm \
    ${refnam} msci noncod.${window} \
    0.5 0.0 1.0 800 \
    bpp_code/popset/${popset}.imap \
    bpp_code/popset/${popset}.pmin \
    ${refnam}/reference/${refnam}.repeats.txt

ls bpp-msci-noncod.${window}/2.vcf2bpp/ | wc -l 


## Run MSCI Models *********************************************************
popset="ilm_msci"; 
cd ${popset}
cd bpp-msci-noncod.${window}
tmodel="modelA"; echo ${tmodel}; mkdir ${tmodel}

#/ 3. Run BPP
sbatch bpp_code/msci/3.msci.RunBPP.slurm \
        bpp_code/ \
        noncod.${window} \
        ${popset} \
        ${tmodel} \
        1 1 1 \
        rep1

# ## Resume MSCi runs that halted 
# popset="mmp_msci"; tmodel="modelA"
# cd $popset/bpp-msci-noncod.${window}
# sbatch bpp_code/msci/4.msci.ResumeBpp.slurm \
#         ${popset} \
#         ${tmodel} \
#         rep6
# cd ../..

