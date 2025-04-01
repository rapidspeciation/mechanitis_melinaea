#!/bin/bash
#SBATCH -n 1                                    # Number of cores requested
#SBATCH -N 1                                    # Ensure that all cores are on one machine
#SBATCH -t 4300                                 # Runtime in minutes
#SBATCH -p sapphire                             # Partition to submit to
#SBATCH --mem-per-cpu=8000                      # Memory per cpu in MB (see also --mem-per-cpu)
#SBATCH --open-mode=append                      #
#SBATCH -o runMSCi_%j.out                       # Standard out goes to this file
#SBATCH -e runMSCi_%j.err                       # Standard err goes to this filehostname



#/ Variables
bppdir=$1

slocus=$2;
popset=$3;
tmodel=$4;

iA=$5;
iB=$6;
iC=$7;

rep=$8;



## Run BPP ******************************************************************************************************
mkdir ${tmodel}
mkdir ${tmodel}/${rep}
cd ${tmodel}/${rep}
cp ../../${slocus}.seqs.txt .


#/ Copy imap file
cp ${bppdir}/popset/${popset}.imap ${popset}.imap.txt

#/ Update control file
cp ${bppdir}/models/msci.${tmodel}.${rep}.ctl msci.${tmodel}.ctl
sed -i "s,seqs,${slocus}.seqs,g" msci.${tmodel}.ctl
sed -i "s,imap,${popset}.imap,g" msci.${tmodel}.ctl
sed -i "s,outp,${popset}.outp,g" msci.${tmodel}.ctl
sed -i "s,mcmc.txt,${popset}.mcmc.txt,g" msci.${tmodel}.ctl

#/ number of loci 
nloci=`ls ../../2.vcf2bpp/* | wc -l`; echo $nloci; export nloci=$nloci
sed -i "s,numberLoci,${nloci},g" msci.${tmodel}.ctl

#/ number of individuals
sed -i "s,ipopA,$iA,g" msci.${tmodel}.ctl
sed -i "s,ipopB,$iB,g" msci.${tmodel}.ctl
sed -i "s,ipopC,$iC,g" msci.${tmodel}.ctl

#/ run bpp 
bpp --cfile msci.${tmodel}.ctl




