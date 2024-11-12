#!/bin/bash
#BSUB -M 25000
#BSUB -R "select[mem>25000] rusage[mem=25000] span[hosts=1]"
#BSUB -n 8
#BSUB -q normal
#BSUB -J align_minimap[7-8]
#BSUB -o log/minimap_lsf_output.%J.%I.txt
#BSUB -e log/minimap_lsf_error.%J.%I.txt

#make a tab sep list of fasta files to align including directories called "list_species.txt"

taxa1=$(sed -n -e ${LSB_JOBINDEX}p list_species.txt |cut -f1)
taxa2=$(sed -n -e ${LSB_JOBINDEX}p list_species.txt |cut -f2)

sp1=$(basename $taxa1)
sp2=$(basename $taxa2)

alignment=${sp1}_${sp2}_${LSB_JOBINDEX}
input_ref=${taxa1}
input_query=${taxa2}
output_paf=${alignment}.paf
output_err=log/${alignment}.err


module load minimap2/2.27--he4a0461_1 

minimap2 \
    -t 8 \
    -x asm10 \
    -o $output_paf \
    $input_ref \
    $input_query \
    2> $output_err
