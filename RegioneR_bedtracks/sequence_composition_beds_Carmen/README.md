### Date created: 24 Sept 2026
### Created by: CP Webster

Scripts and bed files in this directory are:

1. Alignments on Sol

run_minimap2_cpw.sh -> Sbatch script used to generated alignments of ONT-UL and PB-HiFi reads to the assembled Gila Monster T2T genome.

run_mosdepth_cpw.sh -> Sbatch script used to calculate sequence depth in 500kb windows for both ONT-UL and PB-HiFi alignments, which yeilded:

ONT_UL_gila_depth_500k.regions.bed
PB_HiFi_gila_depth_500k.regions.bed

2. Normalize sequence depth

Normalize_SeqDepth.R -> R script to normalize the sequencing depth value for both ONT-UL and PB-HiFi sequencing platforms separately and together. Sequence platform depth values analyzed separately were normalized using depth divided by median values; normalized sequence platform depth values were subtracted ((ONT-UL)-(PB-HiFi)) from one another, which yeilded:

ONT_UL_Normalized_SeqDepth_500kb.bed
PB_HiFi_Normalized_SeqDepth_500kb.bed
NormalizedDifference_SeqDepth_500kb.bed

3. Normalize sequence composition

Normalize_SeqCom.R -> Used Heloderma_suspectum.final.reference.500kb-windows.seqcomp.out to parse GA, GC, and GT values, sum, and normalize by dividing the value by 500000, which yeilded:

GA_Normalized_SeqComp_500kb.bedGC_Normalized_SeqComp_500kb.bedGT_Normalized_SeqComp_500kb.bed

