#!bin/bash

#create bed from text file
tail -n +2 gilamonster_rptmsk_F4.PREDICT.BSCORRECTED.txt | cut -f1-3,5 > gilamonster_rptmsk_F4.PREDICT.BSCORRECTED.bed

#create 500kb windows file, excluding mtDNA
bedtools makewindows -g Heloderma_suspectum.final.reference.fasta.fai -w 500000 | grep -v 'mtDNA' > Heloderma_suspectum.final.reference.500kb_windows.bed

#intersect recombination and 500kb calculating window MEAN
bedtools map -c 4 -o mean -a Heloderma_suspectum.final.reference.500kb_windows.bed -b gilamonster_rptmsk_F4.PREDICT.BSCORRECTED.bed > gilamonster_rptmsk_F4.PREDICT.BSCORRECTED.500kb.bed

#BEWARE, nonexistent windows in query result in "." values, e.g.
chr1    0       500000  .
chr1    500000  1000000 1.52e-11
chr1    1000000 1500000 7.6e-12
chr1    1500000 2000000 0
chr1    2000000 2500000 4.41e-12
chr1    2500000 3000000 2.015e-11
chr1    3000000 3500000 4.465e-12
chr1    3500000 4000000 3.16e-11
chr1    4000000 4500000 1.2145e-10
chr1    4500000 5000000 7.586666667e-11
