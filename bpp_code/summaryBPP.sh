#/ check if files finished or need to be resumed; if done copy final files to finalized folder
mkdir 3.2.bpp/finished/
bash bpp_code/scripts/bpp_status.sh 100001
#/ concatenate output files // tree Counts
rm -f bpp_block_treeCounts.txt
ls 3.2.bpp/finished/*out* | xargs -n 1 -P 1 sh -c '
    perl bpp_code/scripts/parseBPP_output.pl $0 \
    >> bpp_block_treeCounts.txt
'
# transform file to more readable format
sed -i 's,3.2.bpp/finished/,,g' bpp_block_treeCounts.txt
sed -i 's,.out.txt,,g' bpp_block_treeCounts.txt
sed -i 's,\.,\t,g' bpp_block_treeCounts.txt
sed -i 's,block\_,,g' bpp_block_treeCounts.txt
