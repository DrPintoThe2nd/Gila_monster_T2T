#!bin/bash

#create 500kb windows file, excluding mtDNA
#bedtools makewindows -g Heloderma_suspectum.final.reference.fasta.fai -w 500000 | grep -v 'mtDNA' > Heloderma_suspectum.final.reference.500kb_windows.bed

windows="../Heloderma_suspectum.final.reference.500kb_windows.bed"
telos="Heloderma_suspectum.final.reference.telo.bed"
sat_arrays="Heloderma_suspectum.final.reference.aniann.bed"
telo_free="Heloderma_suspectum.final.reference.telo-free_aniann.bed"

bedtools intersect -v -a ${sat_arrays} -b ${telos} > ${telo_free}

for line in telo-free_aniann telo
do
bedtools coverage -mean -a ${windows} -b Heloderma_suspectum.final.reference.$line\.bed > Heloderma_suspectum.final.reference.$line\.500kb.bed
Rscript ../auto_plotting.R Heloderma_suspectum.final.reference.$line\.500kb.bed
done




