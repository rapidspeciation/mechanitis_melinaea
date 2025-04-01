
#/ Arguments
refnam=$1; export refnam=$refnam
region=$2; export region=$region
slocus=$3; export slocus=$slocus
nlocus=$4; export nlocus=$nlocus
breaks=$5; export breaks=$breaks
scalen=$6; export scalen=$scalen

# refnam="melinaea" export refnam=$refnam
# region="noncod"; export region=$region
# slocus="l300_g2000"; export slocus=$slocus
# nlocus=100; export nlocus=$nlocus
# breaks="none"; export breaks=$breaks
# scalen="/n/holyscratch01/mallet_lab/fseixas/ithomiini/$refnam/reference/${refnam}.scaffold.sizes"; export scalen=$scalen

# refnam="mechanitis" export refnam=$refnam
# region="noncod"; export region=$region
# slocus="l300_g2000"; export slocus=$slocus
# nlocus=100; export nlocus=$nlocus
# breaks="none"; export breaks=$breaks
# scalen="/n/holyscratch01/mallet_lab/fseixas/ithomiini/$refnam/reference/${refnam}.scaffold.sizes"; export scalen=$scalen

# refnam="ilMecMess1.1"; export refnam=$refnam
# region="noncod"; export region=$region
# slocus="l300_g2000"; export slocus=$slocus
# nlocus=100; export nlocus=$nlocus
# breaks="none"; export breaks=$breaks
# scalen="/n/home12/fseixas/code/heliconius_seixas/mechanatis-code/ilMecMess1.1.scalen.txt"; export scalen=$scalen


#/ Create directories
mkdir bpp-a01-${region}.${slocus}/3.0.blocks
mkdir bpp-a01-${region}.${slocus}/3.1.blocksAlign
mkdir bpp-a01-${region}.${slocus}/3.2.bpp


## Prepare blocks of loci ******************************************************************************************

#/ List loci to concatenate
cd bpp-a01-${region}.${slocus}/
cd 2.vcf2bpp
ls *.bpp.fasta > ../finalLoci.list
cd ..
sed 's,\.,\t,g' finalLoci.list | sort -k1,1 -k2,2n | sed 's,\t,.,g' > tmp; mv tmp finalLoci.list

#/ list loci by chromosome
for scaff in ` cut -f1 ${scalen} | uniq `; do
    echo $scaff
    egrep -w $scaff finalLoci.list > 3.0.blocks/$scaff.list
done
# rm finalLoci.list

#/ Split by breakpoint (e.g. inversions)
if [[ $breaks = "none" ]]
then
    echo "skip breaking"
else
    echo "splitting chromosomes"
    module load python/3.10.9-fasrc01
    cat ${breaks} | xargs -n 2 sh -c '
        echo $0
        python bpp_code/scripts/breakBlocks.py \
            3.0.blocks/$0.list \
            ${breaks} \
            3.0.blocks/ \
            $0
            rm 3.0.blocks/$0.list 
    '
fi

#/ Determine blocks of loci
for lst in ` ls 3.0.blocks/*.list | sed 's,3.0.blocks/,,g' | sed 's,.list,,g' `; do
    echo $lst
    split -l ${nlocus} -a 3 --numeric-suffixes=1 3.0.blocks/$lst.list 3.0.blocks/$lst.block\_;
    rm -f 3.0.blocks/$lst.list
done
wc -l 3.0.blocks/*


## Collapse or remove blocks at the end of chromosomes when too few loci **************************************
rm -f locusList.txt
ls 3.0.blocks/*  | xargs -n 1 -P 1 sh -c '
    nlocus=`wc -l $0`
    printf "%s\t%s\n" $nlocus >> locusList.txt;
'
awk ' {print $2,$1 }' locusList.txt | sed 's, ,\t,g' > tmp; mv tmp locusList.txt

#/ Note: this code only works when blocks with less than [min] alignments occur at the end of the chromosome (last block) 
for scaff in ` ls 3.0.blocks/ | cut -d'.' -f1 | uniq `; do
    egrep ${scaff} locusList.txt | awk -v min=50 '$2 < min {printf "%s\n%s\n",f,$1} {f=$1}' > toMerge.txt
    # define new name of files to concatenate: 1st in list
    newnam=$(head -n 1 toMerge.txt); echo $newnam
    todrop=$(tail -n 1 toMerge.txt); echo $todrop
    # concatenate files and remove unnecessary ones (last)
    for file in `cat toMerge.txt`; do cat $file >> tmp; done
    mv tmp $newnam
    rm -f $todrop
done
rm -f toMerge.txt

# calculate size of blocks
for file in `ls 3.0.blocks`; do
    nam=` echo $file | sed 's,\.,\t,g' | awk '{print $1}'`
    sta=`head -n 1 3.0.blocks/$file | sed 's/\./\t/g' | awk '{print $2}'`
    end=`tail -n 1 3.0.blocks/$file | sed 's/\./\t/g' | awk '{print $3}'`
    len=`echo $(( ($end-$sta))) `
    printf '%s\t%s\t%s\t%s\n' $nam $sta $end $len
done > blockSummary.txt


#/ Concatenate alignments in the same block ***********************************************************************
ls 3.0.blocks/* | sed 's,3.0.blocks/,,g' | xargs -n 1 -P 4 sh -c '
    echo $0;
    rm -f 3.1.blocksAlign/$0.seqs.txt
    for aln in `cat 3.0.blocks/$0`; do 
        # echo $aln
        cat 2.vcf2bpp/$aln >> 3.1.blocksAlign/$0.seqs.txt;
    done
'



