#!/bin/bash
#This script will produce synteny plots and intermediate chain files from minimap2 paf-files
#to run do: bash run_syntenyPlotter.sh list_of_alignments.txt
#the file list_of_alignments.txt should contain a list with 3 space sep columns path/to/table_file.paf taxa_name1 taxa_name2 

mkdir plots
mkdir intermediate
#mkdir intermediate/ordered/

taxa_pair=$1

while read -r line; do
Rscript ../Syntenyplotter_paf_v4.R ${line}
done < $taxa_pair
