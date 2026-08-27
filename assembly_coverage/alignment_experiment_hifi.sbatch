#!/bin/bash

#SBATCH -c 48                       # number of cores
#SBATCH --mem=250GB                 # specify additional memory to accomodate analysis
#SBATCH -t 1-00:00:00               # wall time (D-HH:MM)
#SBATCH -p public                  # Uses the
#SBATCH -q public                   # Submits job
#SBATCH -o hifi_align.%j.out           # STDOUT (%j = JobId)
#SBATCH -e hifi_align.%j.err           # STDERR (%j = JobId)
#SBATCH --mail-type=ALL             # Send a notification when a job starts, stops, or fails
#SBATCH --mail-user=bjp004@morningside.edu # send-to address
#SBATCH --export=NONE

module purge

#activate conda env
source activate hifiasm

#align 10X assembly and to haploid, PAR-masked reference
mm2plus -t4 -x asm5 -c --eqx --secondary=no Heloderma_suspectum.final.reference.PAR-masked.fasta.gz Heloderma_10X.contigs.sort.1kb.fasta > Heloderma_10X.contigs.sort.1kb.mm2plus.paf  2>/dev/null

#cat both SMRTcells and do math based on average read length
#cat m84082_240313_102849_s2.hifi_reads.fastq.gz m84082_240313_122818_s3.hifi_reads.fastq.gz > hifi_reads.fastq.gz
#1x coverage == 151,875 reads
#seqkit fq2fa -j 40 hifi_reads.fastq.gz -o hifi_reads.fasta.gz

#subset hifi reads using random seeds
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

# define variables
prefix="asm.hifiasm.70x_alignment_final"
hifi_reads='/scratch/bpinto2/Heloderma_2024/hifi_reads_70x.fasta.gz'
HiC1='/scratch/bpinto2/Heloderma_2024/Gila_monster_R1.fastq.gz'
HiC2='/scratch/bpinto2/Heloderma_2024/Gila_monster_R2.fastq.gz'
diploid_ref='/scratch/bpinto2/Heloderma_2024/Heloderma_suspectum.final.reference.diploid.fasta'

hifiasm -o ${prefix} -t48 --h1 ${HiC1} --h2 ${HiC2} ${hifi_reads}

#loop across each subset
for line in 70x 60x 50x 40x 30x 20x
do

#convert gfa to fasta
gfatools gfa2fa asm.hifiasm.$line\_alignment_final.hic.hap1.p_ctg.gfa > asm.hifiasm.$line\_alignment_final.hic.hap1.p_ctg.fa
gfatools gfa2fa asm.hifiasm.$line\_alignment_final.hic.hap2.p_ctg.gfa > asm.hifiasm.$line\_alignment_final.hic.hap2.p_ctg.fa

#merge into diploid assembly
cat asm.hifiasm.$line\_alignment_final.hic.hap1.p_ctg.fa asm.hifiasm.$line\_alignment_final.hic.hap2.p_ctg.fa > asm.hifiasm.$line\_alignment_final.hic.diploid.p_ctg.fa

#align diploid hifi to diploid T2T
mm2plus -t48 -x asm5 -c --eqx --secondary=no ${diploid_ref} asm.hifiasm.$line\_alignment_final.hic.diploid.p_ctg.fa > asm.hifiasm.$line\_alignment_final.hic.diploid.mm2plus.paf  2>/dev/null

#parse out haps 1 & 2 and remove mtDNA
grep "_alt" asm.hifiasm.$line\_alignment_final.hic.diploid.mm2plus.paf > hap2_$line\.paf
grep -v "_alt" asm.hifiasm.$line\_alignment_final.hic.diploid.mm2plus.paf > tmp
grep "chrW" tmp >> hap2_$line\.paf
grep -E -v '16-S_mtDNA|chrW' tmp > hap1_$line\.paf

#check outputs
cut -f6 hap1_$line\.paf | sort | uniq
cut -f6 hap2_$line\.paf | sort | uniq

done

