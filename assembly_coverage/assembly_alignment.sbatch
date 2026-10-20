#!/bin/bash

#SBATCH -c 48                       # number of cores
#SBATCH --mem=100GB                 # specify additional memory to accomodate analysis
#SBATCH -t 0-04:00:00               # wall time (D-HH:MM)
#SBATCH -p htc                  # Uses the
#SBATCH -q public                   # Submits job
#SBATCH -o hifi_align.%j.out           # STDOUT (%j = JobId)
#SBATCH -e hifi_align.%j.err           # STDERR (%j = JobId)
#SBATCH --mail-type=ALL             # Send a notification when a job starts, stops, or fails
#SBATCH --mail-user=bjp004@morningside.edu # send-to address
#SBATCH --export=NONE

module purge

#activate conda env
source activate hifiasm

#cat both SMRTcells and do math based on average read length
#cat m84082_240313_102849_s2.hifi_reads.fastq.gz m84082_240313_122818_s3.hifi_reads.fastq.gz > hifi_reads.fastq.gz
#1x coverage == 151,875 reads
#seqkit fq2fa -j 40 hifi_reads.fastq.gz -o hifi_reads.fasta.gz

######subset hifi reads using random seeds
#60x
#seqtk sample -s187 hifi_reads.fasta.gz 9112500     | bgzip -c > hifi_reads_60x.fasta.gz
#50x
#seqtk sample -s658 hifi_reads_60x.fasta.gz 7593750 | bgzip -c > hifi_reads_50x.fasta.gz
#40x
#seqtk sample -s109 hifi_reads_50x.fasta.gz 6075000 | bgzip -c > hifi_reads_40x.fasta.gz
#30x
#seqtk sample -s10 hifi_reads_40x.fasta.gz 4556250  | bgzip -c > hifi_reads_30x.fasta.gz
#25x
#seqtk sample -s17 hifi_reads_30x.fasta.gz 3796886  | bgzip -c > hifi_reads_25x.fasta.gz
#20x
#seqtk sample -s101 hifi_reads_25x.fasta.gz 3037500 | bgzip -c > hifi_reads_20x.fasta.gz

######define variables
#prefix="asm.hifiasm.70x_alignment_final"
#hifi_reads='/scratch/bpinto2/Heloderma_2024/hifi_reads_70x.fasta.gz'
#HiC1='/scratch/bpinto2/Heloderma_2024/Gila_monster_R1.fastq.gz'
#HiC2='/scratch/bpinto2/Heloderma_2024/Gila_monster_R2.fastq.gz'
#hap1_ref='/scratch/bpinto2/Heloderma_2024/hifi_asm/Heloderma_suspectum_10-3-25.hap1.HiFi-polished.fasta.gz'
#hap2_ref='/scratch/bpinto2/Heloderma_2024/hifi_asm/Heloderma_suspectum_10-3-25.hap2.HiFi-polished.fasta.gz'

######loop across each subset
#for line in 20x 30x 40x 50x 60x
#do
#hifiasm -o ${prefix} -t48 --h1 ${HiC1} --h2 ${HiC2} hifi_reads_$line\.fasta.gz
######convert gfa to fasta
#gfatools gfa2fa asm.hifiasm.$line\_alignment_final.hic.hap1.p_ctg.gfa > asm.hifiasm.$line\_alignment_final.hic.hap1.p_ctg.fa
#gfatools gfa2fa asm.hifiasm.$line\_alignment_final.hic.hap2.p_ctg.gfa > asm.hifiasm.$line\_alignment_final.hic.hap2.p_ctg.fa
######merge into diploid assembly
#cat asm.hifiasm.$line\_alignment_final.hic.hap1.p_ctg.fa asm.hifiasm.$line\_alignment_final.hic.hap2.p_ctg.fa > asm.hifiasm.$line\_alignment_final.hic.diploid.p_ctg.fa
#done

ref='/scratch/bpinto2/Heloderma_2024/hifi_asm/Heloderma_suspectum.final.reference.fasta.gz'

mm2plus -t48 -x asm5 -c --secondary=no --paf-no-hit ${ref} Gila_10X.contigs.fasta.gz > asm.hifiasm.10X_alignment_final.reference.mm2plus.paf

for line in 20x 30x 40x 50x 60x 70x
do
#align diploid hifi assembly to each reference haplotype, then correct haplotype assignments using hiphap
mm2plus -t48 -x asm5 -c --secondary=no --paf-no-hit ${ref} asm.hifiasm.$line\_alignment_final.hic.diploid.p_ctg.fa > asm.hifiasm.$line\_alignment_final.reference.mm2plus.paf
done

source activate scinkd3

for line in 20x 30x 40x 50x 60x 70x
do
#plot the data
Rscript paf_cov_plotting_capped.R asm.hifiasm.$line\_alignment_final.reference.mm2plus.paf

done

do
#align diploid hifi assembly to each reference haplotype, then correct haplotype assignments using hiphap
mm2plus -t48 -x asm5 -c --secondary=no --paf-no-hit ${ref} asm.hifiasm.$line\_alignment_final.hic.diploid.p_ctg.fa > asm.hifiasm.$line\_alignment_final$
done

source activate scinkd3
#plot the HiFi assembly data
for line in 20x 30x 40x 50x 60x 70x
do
for line in 20x 30x 40x 50x 60x 70x
do
Rscript paf_cov_plotting_capped.R asm.hifiasm.$line\_alignment_final.reference.mm2plus.paf
done

#plot the 10X genomics data
Rscript paf_cov_plotting_capped.R asm.hifiasm.10X_alignment_final.reference.mm2plus.paf
