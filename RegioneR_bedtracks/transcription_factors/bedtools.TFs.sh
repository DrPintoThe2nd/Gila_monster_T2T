#!bin/bash

#create 500kb windows file, excluding mtDNA
#bedtools makewindows -g Heloderma_suspectum.final.reference.fasta.fai -w 500000 | grep -v 'mtDNA' > Heloderma_suspectum.final.reference.500kb_windows.bed

windows="Heloderma_suspectum.final.reference.500kb_windows.bed"
TF_in="Hsusp_all_TF_annotations.gff"
TF_cov="Heloderma_suspectum_all_TF_annotations-cov.500kb.bed"
TF_count="Heloderma_suspectum_all_TF_annotations-count.500kb.bed"

#intersect HYPO-methylated
#coverage of sites
bedtools coverage -mean -a ${windows} -b ${TF_in} | grep -v "16-S_mtDNA" > ${TF_cov}
#counts of sites
bedtools coverage -counts -a ${windows} -b ${TF_in} | grep -v "16-S_mtDNA" > ${TF_count}

Rscript ../auto_plotting.R ${TF_cov}
Rscript ../auto_plotting.R ${TF_count}






