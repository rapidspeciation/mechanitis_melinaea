## Run BPP analysis in blocks ******************************************************************************************

## set variables & export
refnam=$1; export refnam=$refnam
region=$2; export region=$region
slocus=$3; export slocus=$slocus

## BPP details
#/ sline: species line (e.g. 6 ismout melgui pargp1 pargp3 elegp3 elegp1)
#/ iline: individuals/pop line (e.g. 1 3 3 3 3 3)
#/ tline1: initial tree rep1 - e.g  "(melgui,(ismout,((pargp1,pargp3),(elegp3,elegp1))));" 
#/ tline2: initial tree rep2 - e.g  "(melgui,(ismout,((pargp1,elegp3),(pargp3,elegp1))));" 
#/ tline3: initial tree rep3 - e.g  "(melgui,(ismout,((pargp1,elegp1),(pargp3,elegp3))));" 
sline=$4; export sline=$sline 
iline=$5; export iline=$iline 
tline=$6; export tline=$tline
phase=$7; export phase=$phase

# sline="8    Hcha_i  Hcha_s  Hper  Hsar  Hsap  Hhor Hera  Hher"
# iline="1  1  1  1  1  1  1  1"
# tline="(((((Hcha_i,Hcha_s),Hper),(Hsar,Hsap)),Hhor),(Hera,Hher));"
# export sline=$sline 
# export iline=$iline 
# export tline1=$tline

imap=$8;  export imap=$imap
name=$9;  export name=$name
repl=${10}; export repl=$repl

echo $repl

## Prepare control files ***********************************************************************************************
for block in `ls 3.0.blocks/* | sed 's,3.0.blocks/,,g'`; do
    echo $block
    mkdir 3.2.bpp/$block.rep${repl}/
    # Copy imap file 
    cp ${imap} 3.2.bpp/$block.rep${repl}/bpp01.imap
    # copy control (.ctl) file 
    cp ~/code/heliconius_seixas/workflows/bpp/bpp01/bpp01.ctl   3.2.bpp/$block.rep${repl}/bpp01.ctl
    # record number of loci per block
    nlocus=`egrep ${name} 3.1.blocksAlign/$block.seqs.txt | wc -l`
    # printf "%s\t%s\n" $block $nlocus >> bpp-a01.loci_summary.txt                                                  # print number of loci
    # change ctl file            
    sed -i "s,numberLoci,$nlocus,g" 3.2.bpp/$block.rep${repl}/bpp01.ctl                                               # specify number of loci
    sed -i "s,seqs.txt,3.1.blocksAlign/$block.seqs.txt,g" 3.2.bpp/$block.rep${repl}/bpp01.ctl                         # specify alignments file
    sed -i "s,imap.txt,3.2.bpp/$block.rep${repl}/bpp01.imap,g" 3.2.bpp/$block.rep${repl}/bpp01.ctl                      # specify imap file
    sed -i "s,outp.txt,3.2.bpp/$block.rep${repl}/$block.rep${repl}.out.txt,g" 3.2.bpp/$block.rep${repl}/bpp01.ctl         # specify output file
    sed -i "s,mcmc.txt,3.2.bpp/$block.rep${repl}/$block.rep${repl}.mcmc.list,g" 3.2.bpp/$block.rep${repl}/bpp01.ctl       # specify mcmc file name
    sed -i "s/sppline/${sline}/g" 3.2.bpp/$block.rep${repl}/bpp01.ctl                                                 # replace species line
    sed -i "s/indline/${iline}/g" 3.2.bpp/$block.rep${repl}/bpp01.ctl                                                 # replace individual count line
    sed -i "s/treline/${tline}/g" 3.2.bpp/$block.rep${repl}/bpp01.ctl                                                 # replace tree line
    sed -i "s/phaseln/${phase}/g" 3.2.bpp/$block.rep${repl}/bpp01.ctl                                                 # replace phase line

    sed -i "s/nburn/50000/g"  3.2.bpp/$block.rep${repl}/bpp01.ctl                                                     # define burnin
    sed -i "s/nfreq/10/g"     3.2.bpp/$block.rep${repl}/bpp01.ctl                                                     # define sampling frequency
    sed -i "s/nstep/100000/g" 3.2.bpp/$block.rep${repl}/bpp01.ctl                                                     # define number of samples
done

