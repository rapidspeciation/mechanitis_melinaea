#!/bin/bash
#BSUB -M 20000
#BSUB -R "select[mem>20000] rusage[mem=20000] span[hosts=1]"
#BSUB -n 12
#BSUB -q long
#BSUB -J busco
#BSUB -o log/busco_lsf_output.%J.txt
#BSUB -e log/busco_lsf_error.%J.txt

#variables
#inputdir for genome fasta, only need to specify directory if running on all genomes
input_dir=../01_data/
output_dir=busco_out

module load busco/5.7.1--pyhdfd78af_0 
#module load busco/5.4.6--pyhdfd78af_0
#usage: busco -i [SEQUENCE_FILE] -l [LINEAGE] -o [OUTPUT_NAME] -m [MODE] [OTHER OPTIONS]
#example for multiple genomes
#busco -i ./protocol3/bact_genomes -l mycoplasmatales_odb10 -m geno -o busco_out_mycoplas_genomes -c 12
#--restart if stopped for some reason

busco -i $input_dir \
    -l lepidoptera_odb10 \
    -m geno \
    -o $output_dir \
    -c 12

