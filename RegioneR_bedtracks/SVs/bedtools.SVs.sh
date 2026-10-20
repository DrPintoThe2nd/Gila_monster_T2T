#!bin/bash

#create 500kb windows file, excluding mtDNA
bedtools makewindows -g Heloderma_suspectum.final.reference.fasta.fai -w 500000 | grep -v 'mtDNA' > Heloderma_suspectum.final.reference.500kb_windows.bed

SVs='Heloderma_suspectum.final.reference.longdust.bed'
windows='Heloderma_suspectum.final.reference.500kb_windows.bed'
out_counts='Heloderma_suspectum.final.reference.longdust.counts.bed'
out_mean='Heloderma_suspectum.final.reference.longdust.mean.bed'

#intersect
#coverage of overlaps
bedtools coverage -mean -a ${windows} -b ${SVs} > ${out_mean}

#counts of predictors
bedtools coverage -counts -a ${windows} -b ${SVs} > ${out_counts}
