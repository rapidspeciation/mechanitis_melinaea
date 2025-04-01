#!/bin/bash
#SBATCH -n 12                                   # Number of cores requested
#SBATCH -N 1                                    # Ensure that all cores are on one machine
#SBATCH -t 600                                  # Runtime in minutes
#SBATCH -p sapphire,shared                      # Partition to submit to
#SBATCH --mem-per-cpu=2000                      # Memory per cpu in MB (see also --mem-per-cpu)
#SBATCH --open-mode=append                      #
#SBATCH -o RM_%j.out                            # Standard out goes to this file
#SBATCH -e RM_%j.err                            # Standard err goes to this filehostname

fasref=$1
export fasref=$fasref

#/ Load libraries
module load python/3.10.12-fasrc01
mamba activate repeatmasker

## run RepeatMasker
RepeatMasker \
    -pa 12 \
    -species Lepidoptera \
    -a \
    -xsmall \
    -nocut \
    -gff \
    ${fasref}
