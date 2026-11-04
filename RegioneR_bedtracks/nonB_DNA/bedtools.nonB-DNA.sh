#!bin/bash

#create 500kb windows file, excluding mtDNA
#bedtools makewindows -g Heloderma_suspectum.final.reference.fasta.fai -w 500000 | grep -v 'mtDNA' > Heloderma_suspectum.final.reference.500kb_windows.bed

windows="../Heloderma_suspectum.final.reference.500kb_windows.bed"

#intersect
for line in ZDNAHunter g4Discovery gfa.APR gfa.IR gfa.TRI
do

#coverage of overlaps
bedtools coverage -mean -a ${windows} -b chrG.Heloderma_suspectum.$line\.bed > chrG.Heloderma_suspectum.$line\.cov.bed

#counts of predictors
bedtools coverage -counts -a ${windows} -b chrG.Heloderma_suspectum.$line\.bed > chrG.Heloderma_suspectum.$line\.counts.bed

Rscript ../auto_plotting.R chrG.Heloderma_suspectum.$line\.counts.bed
Rscript ../auto_plotting.R chrG.Heloderma_suspectum.$line\.cov.bed

done

