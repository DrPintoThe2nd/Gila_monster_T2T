#!bin/bash

#create 500kb windows file, excluding mtDNA
bedtools makewindows -g Heloderma_suspectum.final.reference.fasta.fai -w 500000 | grep -v 'mtDNA' > Heloderma_suspectum.final.reference.500kb_windows.bed

#collect proportion of sites with <0.25 methylation read frequency (considered HYPOmethylated)
awk ' $5 <= 0.25 ' 2_Female_5mC_5hmC_covFiltered4.tsv > 2_Female_5mC_5hmC_covFiltered4.HYPO.bed

#collect proportion of sites with <0.25 methylation read frequency (considered HYPERmethylated)
awk ' $5 >= 0.75 ' 2_Female_5mC_5hmC_covFiltered4.tsv > 2_Female_5mC_5hmC_covFiltered4.HYPER.bed

windows="Heloderma_suspectum.final.reference.500kb_windows.bed"
HYPO_in="2_Female_5mC_5hmC_covFiltered4.HYPO.bed"
HYPER_in="2_Female_5mC_5hmC_covFiltered4.HYPER.bed"
HYPO_cov="Heloderma_suspectum_16-5_5mC_5hmC_covFiltered4.HYPO-cov.500kb.bed"
HYPO_count="Heloderma_suspectum_16-5_5mC_5hmC_covFiltered4.HYPO-count.500kb.bed"
HYPER_cov="Heloderma_suspectum_16-5_5mC_5hmC_covFiltered4.HYPER-cov.500kb.bed"
HYPER_count="Heloderma_suspectum_16-5_5mC_5hmC_covFiltered4.HYPER-count.500kb.bed"

#intersect HYPO-methylated
#coverage of sites
bedtools coverage -mean -a ${windows} -b ${HYPO_in} > ${HYPO_cov}
#counts of sites
bedtools coverage -counts -a ${windows} -b ${HYPO_in} > ${HYPO_count}

#intersect HYPER-methylated
#coverage of sites
bedtools coverage -mean -a ${windows} -b ${HYPER_in} > ${HYPER_cov}
#counts of sites
bedtools coverage -counts -a ${windows} -b ${HYPER_in} > ${HYPER_count}






