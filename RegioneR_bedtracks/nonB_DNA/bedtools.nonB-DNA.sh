#!bin/bash

#create bed from text file
tail -n +2 gilamonster_rptmsk_F4.PREDICT.BSCORRECTED.txt | cut -f1-3,5 > gilamonster_rptmsk_F4.PREDICT.BSCORRECTED.bed

#create 500kb windows file, excluding mtDNA
bedtools makewindows -g Heloderma_suspectum.final.reference.fasta.fai -w 500000 | grep -v 'mtDNA' > Heloderma_suspectum.final.reference.500kb_windows.bed

#intersect
for line in ZDNAHunter g4Discovery gfa.APR gfa.IR gfa.TRI
do

#coverage of overlaps
bedtools coverage -mean -a Heloderma_suspectum.final.reference.500kb_windows.bed -b chrG.Heloderma_suspectum.$line\.bed > chrG.Heloderma_suspectum.$line\.mean.bed

#counts of predictors
bedtools coverage -counts -a Heloderma_suspectum.final.reference.500kb_windows.bed -b chrG.Heloderma_suspectum.$line\.bed > chrG.Heloderma_suspectum.$line\.counts.bed

done
