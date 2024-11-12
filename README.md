# Code and data underlying the manuscript by van der Heijden et al. on _Mechanitis_ and _Melinaea_ butterfly diversification.

https://www.biorxiv.org/content/10.1101/2024.07.07.602206v1

**Processing of raw reads**
- [Snakemake file 1](https://github.com/rapidspeciation/mechanitis_melinaea/blob/main/scripts/snakefile%20step1)
- [Snakemake file 2](https://github.com/rapidspeciation/mechanitis_melinaea/blob/main/scripts/snakefile%20step2)

**Figure 1:**
- [Nuclear and mitochondrial phylogenetic tree](https://github.com/rapidspeciation/mechanitis_melinaea/blob/main/scripts/phylogenetics)
- [R script to plot the co-phylogeny](https://github.com/rapidspeciation/mechanitis_melinaea/blob/main/scripts/cophylo%20plot.R)

**Figure 2:**
- [Calibrated tree](https://github.com/rapidspeciation/mechanitis_melinaea/blob/main/scripts/phylogenetics)
- [R script to plot the distribution maps](https://github.com/rapidspeciation/mechanitis_melinaea/blob/main/scripts/distribution_maps.R)
- [DSuite analysis](https://github.com/rapidspeciation/mechanitis_melinaea/blob/main/scripts/DSuite.sh)
- [TWISST and WindowTree concordance factor](https://github.com/rapidspeciation/mechanitis_melinaea/blob/main/scripts/WindowTree%2BTWISST.sh)

**Figure 3:**
- [Bash script to generate the genome scans with PBS and fd](https://github.com/rapidspeciation/mechanitis_melinaea/tree/main/scripts/nesaea_introgression_PBS.sh)
- [R script to plot the genome scans](https://github.com/rapidspeciation/mechanitis_melinaea/tree/main/Rscript_Mech.nesaea.r)
- Underlying data: [Fst and Dxy values](https://github.com/rapidspeciation/mechanitis_melinaea/tree/main/input/Mechanitis.nesaea.filtered.2000.Fst.Dxy.pi.csv), [introgression values](https://github.com/rapidspeciation/mechanitis_melinaea/tree/main/input/Mechanitis.nesaea.filtered.2000.fd.polBr_nes_lysBr_mess.csv) and [PBS values](https://github.com/rapidspeciation/mechanitis_melinaea/tree/main/input/polW_polB_nes.pbs)

**Figure 4:**
- [R script to filter BUSCO output](https://github.com/rapidspeciation/mechanitis_melinaea/tree/main/scripts/filter_busco_output.R)
- [R script to calculate and plot rearrangements](https://github.com/rapidspeciation/mechanitis_melinaea/tree/main/scripts/rearrangement_stats.Rmd)
- [R script to detect breakpoints in Mechanitis and Melinaea](https://github.com/rapidspeciation/mechanitis_melinaea/tree/main/scripts/breakpoints.Rmd)
- [R script to plot ribbon plot with whole genome alignments](https://github.com/rapidspeciation/mechanitis_melinaea/tree/main/scripts/Syntenyplotter_paf_multi.R)
- [R script to refine ribbon plot with whole genome alignments](https://github.com/rapidspeciation/mechanitis_melinaea/tree/main/scripts/Syntenyplotter_paf_refining_multi.R)

**Supplementary figures**
- [R script for preliminary synteny ribbon plot with busco genes](https://github.com/rapidspeciation/mechanitis_melinaea/tree/main/scripts/plot_synteny.R)
- [R script for refined synteny ribbon plot with busco genes](https://github.com/rapidspeciation/mechanitis_melinaea/tree/main/scripts/plot_synteny_fine_tuning.R)
- [R script to plot genomescans with breakpoints and boxplots Mech](https://github.com/rapidspeciation/mechanitis_melinaea/tree/main/scripts/genomescan_mec.Rmd)
- [R script to plot genomescans with breakpoints and boxplots Mel](https://github.com/rapidspeciation/mechanitis_melinaea/tree/main/scripts/genomescan_mel.Rmd)
- [bash wrapper script to plot ribbon plot with whole genome alignments between haplotypes](https://github.com/rapidspeciation/mechanitis_melinaea/tree/main/scripts/run_syntenyplotter.R)
- [R script to plot ribbon plot with whole genome alignments between haplotypes](https://github.com/rapidspeciation/mechanitis_melinaea/tree/main/scripts/Syntenyplotter_paf_v4.R)

