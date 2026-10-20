library(tidyverse)

# read in genome wide sequence composition
All_SeqComp <- read.table("Heloderma_suspectum.final.reference.500kb-windows.seQcomp.out", 
                          header = TRUE, sep = "\t")

# name columns
colnames(All_SeqComp) <- c('SC_chr', 'SC_start', 'SC_stop', 'AT', 'GC', 'A', 'C', 
                           'G', 'T', 'N', 'other', 'len_explored')

### GA ###

# parse G and A columns
GA <- All_SeqComp[c('SC_chr', 'SC_start', 'SC_stop', 'A', 'G')]

# add G and A occurrences together 
GA$sum_GA <- GA$A + GA$G

# NORMALIZE -> divide GA sum by 500,000
GA <- GA %>% 
  mutate(sum_GA_norm = sum_GA/500000)

# write to bed
write.table(GA, "GA_Normalized_SeqComp_500kb.bed", 
            sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

### GC ###

# parse G and C columns
GC <- All_SeqComp[c('SC_chr', 'SC_start', 'SC_stop', 'C', 'G')]

# add G and C occurrences together 
GC$sum_GC <- GC$C + GC$G

# NORMALIZE -> divide GC sum by 500,000
GC <- GC %>% 
  mutate(sum_GC_norm = sum_GC/500000)

# write to bed
write.table(GC, "GC_Normalized_SeqComp_500kb.bed", 
            sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

### GT ###

# parse G and T columns
GT <- All_SeqComp[c('SC_chr', 'SC_start', 'SC_stop', 'G', 'T')]

# add G and T occurrences together 
GT$sum_GT <- GT$T + GT$G

# NORMALIZE -> divide GT sum by 500,000
GT <- GT %>% 
  mutate(sum_GT_norm = sum_GT/500000)

# write to bed
write.table(GT, "GT_Normalized_SeqComp_500kb.bed", 
            sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
