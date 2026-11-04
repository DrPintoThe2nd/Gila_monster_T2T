#!bin/bash

#create 500kb windows file, excluding mtDNA
bedtools makewindows -g Heloderma_suspectum.final.reference.fasta.fai -w 500000 | grep -v 'mtDNA' > Heloderma_suspectum.final.reference.500kb_windows.bed

pigz -dc hSus_resolved_all.bed.gz | cut -f7 | sort | uniq > hSus_resolved_types.txt
pigz -dc hSus_resolved_all.bed.gz | awk ' $7 == "DNA" ' > hSus_resolved_DNA.bed
pigz -dc hSus_resolved_all.bed.gz | awk ' $7 == "LINE" ' > hSus_resolved_LINE.bed
pigz -dc hSus_resolved_all.bed.gz | awk ' $7 == "LTR" ' > hSus_resolved_LTR.bed
pigz -dc hSus_resolved_all.bed.gz | awk ' $7 == "RC" ' > hSus_resolved_RC.bed
pigz -dc hSus_resolved_all.bed.gz | awk ' $7 == "SINE" ' > hSus_resolved_SINE.bed

windows="../Heloderma_suspectum.final.reference.500kb_windows.bed"

for line in DNA LINE LTR RC SINE
do
bedtools coverage -mean -a ${windows} -b hSus_resolved_$line\.bed > hSus_resolved_$line\.500kb.bed
Rscript ../auto_plotting.R hSus_resolved_$line\.500kb.bed
done




