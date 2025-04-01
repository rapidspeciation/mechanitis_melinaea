#!/bin/bash
#SBATCH -n 12                                   # Number of cores requested
#SBATCH -N 1                                    # Ensure that all cores are on one machine
#SBATCH -t 300                                  # Runtime in minutes
#SBATCH -p sapphire,shared                      # Partition to submit to
#SBATCH --mem-per-cpu=2000                      # Memory per cpu in MB (see also --mem-per-cpu)
#SBATCH --open-mode=append                      #
#SBATCH -o infoDen_%j.out                       # Standard out goes to this file
#SBATCH -e infoDen_%j.err                       # Standard err goes to this filehostname


#/ Load libraries
module load python/3.10.12-fasrc01


## Arguments ********************************************************************************************************
# regions parameters
refnam=$1;        export refnam=$refnam
buffer=$2;        export buffer=$buffer
wsize=$3;         export wsize=$wsize
wstep=$4;         export wstep=$wstep
lsize=$5;         export lsize=$lsize
lstep=$6;         export lstep=$lstep

vcfdir=$7;         export vcfdir=$vcfdir
ilist=$8;          export ilist=$ilist
cat ${ilist}

#/ Reference genome coordinates
exons="../reference/${refnam}.exons.bed"
scalen="../reference/${refnam}.scaffold.sizes"

#/ Create directories
mkdir bedFiles
mkdir bedFiles/infoDensity/
mkdir 0.vcfflt/



## 1. Genomic coordinates ************************************************************************************

#/ Copy exons bed file
cp $exons bedFiles/${refnam}.exonic.bed

#/ Create buffer around exons: for non coding regions
bedtools slop \
    -b ${buffer} \
    -i ${exons} \
    -g ${scalen} \
    | sort -k1,1V -k2,2n > bedFiles/${refnam}.${buffer}bp_exonbuff.sorted.bed
bedtools merge \
    -i bedFiles/${refnam}.${buffer}bp_exonbuff.sorted.bed \
    > bedFiles/${refnam}.${buffer}bp_exonbuff.merged.bed

# regions outside the buffered exons
bedtools complement \
    -i bedFiles/${refnam}.${buffer}bp_exonbuff.merged.bed \
    -g ${scalen} \
    > bedFiles/${refnam}.noncod.bed
# remove intermediary files
rm -f bedFiles/${refnam}.${buffer}bp_exonbuff.sorted.bed
rm -f bedFiles/${refnam}.${buffer}bp_exonbuff.merged.bed


## check stage
#/ exons and noncod should not overlap
bedtools intersect -a bedFiles/${refnam}.noncod.bed -b $exons
#/ length of noncod
awk '{sum=sum+$3-$2} END {print sum}' bedFiles/${refnam}.noncod.bed




## 2. Define regions to analyze *******************************************************************************

## Subset to individuals and filter out missing genotypes -----
cut -f1 ${scalen} | xargs -n 1 -P 6 sh -c '
    echo $0; 
    # subset
    bcftools view --samples-file ${ilist} -o 0.vcfflt/$0.vcf.gz -O z --threads 2 ${vcfdir}/$0.vcf.gz
    # index
    tabix -p vcf 0.vcfflt/$0.vcf.gz
'


## Regions of low/high density information -----------
#/ Get information density in sliding windows
cat ${scalen} | xargs -n 2 -P 12 sh -c '
    echo $0;
    python bpp_code/scripts/callableRegions.py \
        --vcf_file 0.vcfflt/$0.vcf.gz \
        --chromosome $0 --chrom_size $1 \
        --window-size ${wsize} --step-size ${wstep} \
        > bedFiles/infoDensity/$0.subsample.w${wsize}s${wstep}.txt
'
# high
rm -f bedFiles/wgenome.subsample.w${wsize}s${wstep}.flt025_hig.txt
for scaffold in `cut -f1 ${scalen} `; do
    echo $scaffold
    tail -n+2 bedFiles/infoDensity/${scaffold}.subsample.w${wsize}s${wstep}.txt | \
    awk ' {OFS = "\t" } { if ($4 >= 0.25) print $1,$2,$3 } ' \
    >> bedFiles/wgenome.subsample.w${wsize}s${wstep}.flt025_hig.txt
done
# low
rm -f bedFiles/wgenome.subsample.w${wsize}s${wstep}.flt025_low.txt
for scaffold in `cut -f1 ${scalen} `; do
    echo $scaffold
    tail -n+2 bedFiles/infoDensity/${scaffold}.subsample.w${wsize}s${wstep}.txt | \
    awk ' {OFS = "\t" } { if ($4 < 0.25) print $1,$2,$3 } ' \
    > bedFiles/wgenome.subsample.w${wsize}s${wstep}.flt025_low.txt
done
# intersect low and high
bedtools intersect \
    -a bedFiles/wgenome.subsample.w${wsize}s${wstep}.flt025_hig.txt \
    -b bedFiles/wgenome.subsample.w${wsize}s${wstep}.flt025_low.txt \
    > bedFiles/wgenome.subsample.w${wsize}s${wstep}.flt025_higlow.txt
# remove intersection from high
bedtools subtract \
    -a bedFiles/wgenome.subsample.w${wsize}s${wstep}.flt025_hig.txt \
    -b bedFiles/wgenome.subsample.w${wsize}s${wstep}.flt025_higlow.txt |
    bedtools merge > bedFiles/wgenome.subsample.w${wsize}s${wstep}.flt025_higflt.txt

#/ get only high density information regions 
bedtools intersect \
    -a bedFiles/${refnam}.noncod.bed \
    -b bedFiles/wgenome.subsample.w${wsize}s${wstep}.flt025_higflt.txt |
    sort -k1,1 -k2,2n \
    > bedFiles/wgenome.noncod.callable.bed
bedtools intersect \
    -a bedFiles/${refnam}.exonic.bed \
    -b bedFiles/wgenome.subsample.w${wsize}s${wstep}.flt025_higflt.txt |
    sort -k1,1 -k2,2n \
    > bedFiles/wgenome.exonic.callable.bed

#/ make windows
python bpp_code/scripts/makewindows.py \
    bedFiles/wgenome.noncod.callable.bed \
    bedFiles/wgenome.noncod.l${lsize}_g${lstep}.bed \
    --window_size ${lsize} \
    --window_step ${lstep}
python bpp_code/scripts/makewindows.py \
    bedFiles/wgenome.exonic.callable.bed \
    bedFiles/wgenome.exonic.l${lsize}_g${lstep}.bed \
    --window_size ${lsize} \
    --window_step ${lstep}
