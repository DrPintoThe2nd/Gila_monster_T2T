#!bin/bash

#create 500kb windows file, excluding mtDNA
#bedtools makewindows -g Heloderma_suspectum.final.reference.fasta.fai -w 500000 | grep -v 'mtDNA' > Heloderma_suspectum.final.reference.500kb_windows.bed

windows="Heloderma_suspectum.final.reference.500kb_windows.bed"
gff_in="Heloderma_suspectum.final.reference.gff3"
gff_gene="Heloderma_suspectum.final.reference_gene.gff"
gff_bed="Heloderma_suspectum.final.reference.bed"
gff_cov="Heloderma_suspectum.final.reference-cov.500kb.bed"
gff_count="Heloderma_suspectum.final.reference-count.500kb.bed"

#subset genes from sub-features
awk ' $3 == "gene" ' ${gff_in} > ${gff_gene}

#convert to bed file
gff2bed < ${gff_gene} | grep -v "16-S_mtDNA" > ${gff_bed}

#intersect 
#coverage of sites
bedtools coverage -mean -a ${windows} -b ${gff_bed} | grep -v "16-S_mtDNA" > ${gff_cov}
#counts of sites
bedtools coverage -counts -a ${windows} -b ${gff_bed} | grep -v "16-S_mtDNA" > ${gff_count}

Rscript ../auto_plotting.R ${gff_cov}
Rscript ../auto_plotting.R ${gff_count}






