library(tidyverse)

# read in sequence depth statistics from mosdepth output
pb_depth <- read.table("PB_HiFi_gila_depth_500k.regions.bed", header = FALSE, sep = "\t")
ont_depth <- read.table("ONT_UL_gila_depth_500k.regions.bed", header = FALSE, sep = "\t")

# add column names
colnames(pb_depth)<- c("chr","start","stop", "pb_depth")
colnames(ont_depth)<- c("chr","start","stop", "ont_depth")

# -----------------------------------------------------------------
# WRITE dataframes to bed with normalized sequence depth column
# -----------------------------------------------------------------

# median for PB
pb_depth <- pb_depth %>% 
  mutate(median_depth = median(pb_depth))

# median for ONT
ont_depth <- ont_depth %>% 
  mutate(median_depth = median(ont_depth))

# seq_depth/median_depth PB
pb_depth <- pb_depth %>% 
  mutate(pb_depth_norm = pb_depth/median_depth)

# seq_depth/median_depth ONT
ont_depth <- ont_depth %>% 
  mutate(ont_depth_norm = ont_depth/median_depth)

# write to bed

# ONT-UL
write.table(ont_depth, "ONT_UL_Normalized_SeqDepth_500kb.bed", 
            sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

# PB-HiFi
write.table(pb_depth, "PB_HiFi_Normalized_SeqDepth_500kb.bed", 
            sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

# -----------------------------------------------------------------
# Compare sequence platforms: (ONT-UL norm)-(PB-HiFi norm)
# -----------------------------------------------------------------

# parse needed columns and merge dataframes into one
pb_depth <- pb_depth[c(6)]
ont_depth <- ont_depth[c(1, 2, 3, 6)]
all_depth <- cbind(ont_depth, pb_depth)

# add column with ONT-PB normalized difference
all_depth$difference <- all_depth$ont_depth_norm - all_depth$pb_depth_norm

# write to bed
write.table(all_depth, "NormalizedDifference_SeqDepth_500kb.bed", 
            sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

