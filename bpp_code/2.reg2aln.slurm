#!/bin/bash
#SBATCH -n 16                                    # Number of cores requested
#SBATCH -N 1                                     # Ensure that all cores are on one machine
#SBATCH -t 300                                   # Runtime in minutes
#SBATCH -p sapphire,shared                       # Partition to submit to
#SBATCH --mem-per-cpu=1000                       # Memory per cpu in MB (see also --mem-per-cpu)
#SBATCH --open-mode=append                       #
#SBATCH -o reg2aln_%j.out                        # Standard out goes to this file
#SBATCH -e reg2aln_%j.err                        # Standard err goes to this filehostname

#/ load modules
module load python/3.10.12-fasrc01

#/ Arguments
refnam=$1;    export refnam=$refnam
bpptyp=$2;    export bpptyp=$bpptyp
slocus=$3;    export slocus=$slocus

maxHor=$4;    export maxHor=$maxHor
maxVer=$5;    export maxVer=$maxVer
maxIUP=$6;    export maxIUP=$maxIUP
minPos=$7;    export minPos=$minPos

indmap=$8;    export indmap=$indmap
popmin=$9;    export popmin=$popmin
tomask=${10}; export tomask=$tomask 

#/ Create directories
mkdir bpp-${bpptyp}-${slocus}
mkdir bpp-${bpptyp}-${slocus}/1.vcfRegions
mkdir bpp-${bpptyp}-${slocus}/2.vcf2bpp


#/ Extract, simplify and convert regions // from vcf to fasta
cat bedFiles/wgenome.${slocus}.bed | xargs -n 3 -P 16 sh -c '

    #/ extract region
    bcftools view \
        -r $0:$1-$2 \
        --targets-file ^$tomask \
        0.vcfflt/$0.vcf.gz \
        --output-type z \
        --output-file bpp-${bpptyp}-${slocus}/1.vcfRegions/$0.$1.$2.vcf.gz

    #/ simplify
    plink \
        --vcf bpp-${bpptyp}-${slocus}/1.vcfRegions/$0.$1.$2.vcf.gz \
        --allow-extra-chr \
        --recode vcf-iid \
        --out bpp-${bpptyp}-${slocus}/1.vcfRegions/$0.$1.$2.simple
    rm -f bpp-${bpptyp}-${slocus}/1.vcfRegions/$0.$1.$2.vcf.gz
    rm -f bpp-${bpptyp}-${slocus}/1.vcfRegions/$0.$1.$2.simple.*log
    rm -f bpp-${bpptyp}-${slocus}/1.vcfRegions/$0.$1.$2.simple.*sex

    #/ vcf2fasta
    perl bpp_code/scripts/vcf2fas.iupac.pl \
        bpp-${bpptyp}-${slocus}/1.vcfRegions/$0.$1.$2.simple.vcf \
        bpp-${bpptyp}-${slocus}/2.vcf2bpp/$0.$1.$2.fasta

    #/ filter alignments
    python bpp_code/scripts/fasta_filter.py \
        bpp-${bpptyp}-${slocus}/2.vcf2bpp/$0.$1.$2.fasta \
        bpp-${bpptyp}-${slocus}/2.vcf2bpp/$0.$1.$2.filter.fasta \
        ${maxHor} ${maxVer} ${maxIUP} \
        ${indmap} \
        ${popmin}

    #/ fasta to bpp format
    perl bpp_code/scripts/fasta2bpp.pl \
        bpp-${bpptyp}-${slocus}/2.vcf2bpp/$0.$1.$2.filter.fasta \
        bpp-${bpptyp}-${slocus}/2.vcf2bpp/$0.$1.$2.bpp.fasta \
        $0.$1.$2 \
        ${minPos}
    rm -f bpp-${bpptyp}-${slocus}/2.vcf2bpp/$0.$1.$2.fasta
    rm -f bpp-${bpptyp}-${slocus}/2.vcf2bpp/$0.$1.$2.filter.fasta
'

rm -f bpp-${bpptyp}-${slocus}/1.vcfRegions/*temp*


## concatenate alignments
rm -f bpp-${bpptyp}-${slocus}/${slocus}.seqs.txt
touch bpp-${bpptyp}-${slocus}/${slocus}.seqs.txt
ls bpp-${bpptyp}-${slocus}/2.vcf2bpp/*bpp.fasta | xargs -n 1 sh -c '
    cat $0 >> bpp-${bpptyp}-${slocus}/${slocus}.seqs.txt
'
