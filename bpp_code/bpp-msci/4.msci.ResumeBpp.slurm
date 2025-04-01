#!/bin/bash
#SBATCH -n 1                                    # Number of cores requested
#SBATCH -N 1                                    # Ensure that all cores are on one machine
#SBATCH -t 4300                                 # Runtime in minutes
#SBATCH -p sapphire                             # Partition to submit to
#SBATCH --mem-per-cpu=8000                      # Memory per cpu in MB (see also --mem-per-cpu)
#SBATCH --open-mode=append                      #
#SBATCH -o msciResume_%j.out                    # Standard out goes to this file
#SBATCH -e msciResume_%j.err                    # Standard err goes to this filehostname


# input variables
popset=$1
tmodel=$2
replNb=$3

cd ${tmodel}/${replNb}/

# get last checkpoint file
checkp=`ls ${popset}*.chk | sort -V | tail -n 1`;
echo $checkp

# resume run
bpp --resume $checkp



