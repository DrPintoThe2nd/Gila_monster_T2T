#!bin/bash

#create 500kb windows file, excluding mtDNA
#bedtools makewindows -g Heloderma_suspectum.final.reference.fasta.fai -w 500000 | grep -v 'mtDNA' > Heloderma_suspectum.final.reference.500kb_windows.bed

windows="Heloderma_suspectum.final.reference.500kb_windows.bed"
TF_in="Hsusp_all_TF_annotations.gff"
TF_gene="Hsusp_all_TF_annotations_gene.gff"
TF_bed="Hsusp_all_TF_annotations_gene.bed"
TF_cov="Heloderma_suspectum_all_TF_annotations-cov.500kb.bed"
TF_count="Heloderma_suspectum_all_TF_annotations-count.500kb.bed"

#subset genes from sub-features
awk ' $3 == "gene" ' ${TF_in} > ${TF_gene}

#convert to bed file
gff2bed < ${TF_gene} > ${TF_bed}

#intersect 
#coverage of sites
bedtools coverage -mean -a ${windows} -b ${TF_bed} | grep -v "16-S_mtDNA" > ${TF_cov}
#counts of sites
bedtools coverage -counts -a ${windows} -b ${TF_bed} | grep -v "16-S_mtDNA" > ${TF_count}

Rscript ../auto_plotting.R ${TF_cov}
Rscript ../auto_plotting.R ${TF_count}






