#!/bin/bash

#SBATCH -c 16                  		# number of cores
#SBATCH --mem=48GB              	# specify additional memory to accomodate bbmap
#SBATCH -t 0-02:00:00          	    # wall time (D-HH:MM)
#SBATCH -p htc      	        # Uses the
#SBATCH -q public     	            # Submits job
#SBATCH -o salmon.%j.out       	 	# STDOUT (%j = JobId)
#SBATCH -e salmon.%j.err       	 	# STDERR (%j = JobId)
#SBATCH --mail-type=ALL             # Send a notification when a job starts, stops, or fails
#SBATCH --mail-user=bjp004@morningside.edu # send-to address
#SBATCH --export=NONE

module purge

source activate salmon

#check potential duplication rates
#ls *_1.fastq.gz > samples.txt
#sed -i 's/_1.fastq.gz//' samples.txt
#cat samples.txt | while read line
#do
#echo "dedupe *_[1/2]_val_[1/2].fq.gz..."
#clumpify.sh -Xmx36g in=$line\_1.fastq.gz in2=$line\_2.fastq.gz out=$line\_R1.fastq.gz out2=$line\_R2.fastq.gz dedupe
#done

#fastqc -t 12 *.fastq.gz
#multiqc .

#use longdust LCR annotations for genome soft-masking for Heloderma and Gallus
#Gallus genome missing the PAR, versions: FULL vs. chrW-masked
#longdust GCF_016699485.2_bGalGal1.mat.broiler.GRCg7b_genomic.fasta > GCF_016699485.2_bGalGal1.mat.broiler.GRCg7b_genomic.longdust.bed
#bedtools maskfasta -fi GCF_016699485.2_bGalGal1.mat.broiler.GRCg7b_genomic.fasta -bed GCF_016699485.2_bGalGal1.mat.broiler.GRCg7b_genomic.longdust.bed -fo GCF_016699485.2_bGalGal1.mat.broiler.GRCg7b_genomic.masked.fasta
#bedtools maskfasta -fi GCF_016699485.2_bGalGal1.mat.broiler.GRCg7b_genomic.masked.fasta -bed GCF_016699485.2_bGalGal1.mat.broiler.GRCg7b_genomic.masked.chrW.bed -fo GCF_016699485.2_bGalGal1.mat.broiler.GRCg7b_genomic.masked.chrW_hard-masked.fasta

#make transcripts for each reference
#gffread -g GCF_016699485.2_bGalGal1.mat.broiler.GRCg7b_genomic.fasta -x GCF_016699485.2_bGalGal1.mat.broiler.GRCg7b_transcripts.fasta --no-pseudo GCF_016699485.2_bGalGal1.mat.broiler.GRCg7b_genomic.gff3
#cat GCF_016699485.2_bGalGal1.mat.broiler.GRCg7b_transcripts.fasta GCF_016699485.2_bGalGal1.mat.broiler.GRCg7b_genomic.masked.fasta > GCF_016699485.2_bGalGal1.mat.broiler.GRCg7b_transcripts.salmon-index-plus-decoys.fasta
#gffread -g GCF_016699485.2_bGalGal1.mat.broiler.GRCg7b_genomic.masked.chrW_hard-masked.fasta -x GCF_016699485.2_bGalGal1.mat.broiler.GRCg7b_transcripts_chrW_hard-masked.fasta --no-pseudo GCF_016699485.2_bGalGal1.mat.broiler.GRCg7b_genomic_chrW_hard-masked.gff3
#cat GCF_016699485.2_bGalGal1.mat.broiler.GRCg7b_transcripts_chrW_hard-masked.fasta GCF_016699485.2_bGalGal1.mat.broiler.GRCg7b_genomic.masked.chrW_hard-masked.fasta > GCF_016699485.2_bGalGal1.mat.broiler.GRCg7b_transcripts_chrW_hard-masked.salmon-index-plus-decoys.fasta

#make index for each set of transcripts
#female index (FULL, PAR missing)
#salmon index -p 12 -k 31 -t GCF_016699485.2_bGalGal1.mat.broiler.GRCg7b_transcripts.salmon-index-plus-decoys.fasta -i GRCg7b_transcripts.salmon-index-plus-decoys.idx -d Gallus_gallus_chrnames.txt
#male index (W masked)
#salmon index -p 12 -k 31 -t GCF_016699485.2_bGalGal1.mat.broiler.GRCg7b_transcripts_chrW_hard-masked.salmon-index-plus-decoys.fasta -i GRCg7b_transcripts_chrW_hard-masked.salmon-index-plus-decoys.idx -d Gallus_gallus_chrnames.txt

#longdust run using default params for consortium
#bedtools maskfasta -fi Heloderma_suspectum.final.reference.chrW-masked.fasta -bed Heloderma_suspectum.final.reference.longdust.bed -fo Heloderma_suspectum.final.reference.masked.chrW_hard-masked.fasta 
#bedtools maskfasta -fi Heloderma_suspectum.final.reference.PAR-masked.fasta -bed Heloderma_suspectum.final.reference.longdust.bed -fo Heloderma_suspectum.final.reference.masked.PAR_hard-masked.fasta

#grab transcripts after fixing PAR-masked GFF file (chrW gene IDs are replicates of autosomal gene IDs)
#gffread -g Heloderma_suspectum.final.reference.PAR-masked.fasta -x Heloderma_suspectum.final.reference.PAR-masked.transcripts.fasta --no-pseudo Heloderma_suspectum.final.reference.PAR-masked.fixed.gff3

#construct salmon transcript index inputs (transcriptome + genome/decoys)
#cat Heloderma_suspectum.final.reference.chrW-masked.transcripts.fasta Heloderma_suspectum.final.reference.masked.chrW_hard-masked.fasta > Heloderma_suspectum.final.reference.chrW-masked.salmon-index-plus-decoys.fasta
#cat Heloderma_suspectum.final.reference.PAR-masked.transcripts.fasta Heloderma_suspectum.final.reference.masked.PAR_hard-masked.fasta > Heloderma_suspectum.final.reference.PAR-masked.salmon-index-plus-decoys.fasta

#female index (PAR-masked)
#salmon index -p 12 -k 31 -t Heloderma_suspectum.final.reference.PAR-masked.salmon-index-plus-decoys.fasta -i Heloderma_suspectum.final.reference.PAR-masked.salmon-index-plus-decoys.idx -d Heloderma_suspectum.chrnames.txt
#male index (W masked)
#salmon index -p 12 -k 31 -t Heloderma_suspectum.final.reference.chrW-masked.salmon-index-plus-decoys.fasta -i Heloderma_suspectum.final.reference.chrW-masked.salmon-index-plus-decoys.idx -d Heloderma_suspectum.chrnames.txt

#calc salmon Gila males
#ls SRR*_Gila_*M_R1.fastq.gz > males.txt
#sed -i 's/_R1.fastq.gz//g' males.txt
#sed -i 's/@//' males.txt
#cat males.txt | while read line
#do
#echo "doing salmon stuff..."
#salmon quant -p 16 -i Heloderma_suspectum.final.reference.chrW-masked.salmon-index-plus-decoys.idx -l ISF -1 $line\_R1.fastq.gz -2 $line\_R2.fastq.gz --softclip --writeMappings -g Heloderma_suspectum.final.reference.chrW-masked.gff3 -o $line\.quant | samtools sort -@16 -O BAM -o $line\.bam
#done

#calc salmon Gila females
#ls SRR*_Gila_*F_R1.fastq.gz > females.txt
#sed -i 's/_R1.fastq.gz//g' females.txt
#sed -i 's/@//' females.txt
#cat females.txt | while read line
#do
#echo "doing salmon stuff..."
#salmon quant -p 16 -i Heloderma_suspectum.final.reference.PAR-masked.salmon-index-plus-decoys.idx -l ISF -1 $line\_R1.fastq.gz -2 $line\_R2.fastq.gz --softclip --writeMappings -g Heloderma_suspectum.final.reference.PAR-masked.fixed.gff3 -o $line\.quant | samtools sort -@16 -O BAM -o $line\.bam
#done

#calc salmon Gallus males
ls SRR*_Gallus-male_R1.fastq.gz > Gallus_males.txt
sed -i 's/_R1.fastq.gz//g' Gallus_males.txt
sed -i 's/@//' Gallus_males.txt
cat Gallus_males.txt | while read line
do
echo "doing salmon stuff..."
salmon quant -p 16 -i GRCg7b_transcripts_chrW_hard-masked.salmon-index-plus-decoys.idx -l ISF -1 $line\_R1.fastq.gz -2 $line\_R2.fastq.gz --softclip --writeMappings -g GCF_016699485.2_bGalGal1.mat.broiler.GRCg7b_genomic_chrW_hard-masked.gff3 -o $line\.quant | samtools sort -@16 -O BAM -o $line\.bam
done

#calc salmon Gallus females
ls SRR*_Gallus-female_R1.fastq.gz > Gallus_females.txt
sed -i 's/_R1.fastq.gz//g' Gallus_females.txt
sed -i 's/@//' Gallus_females.txt
cat Gallus_females.txt | while read line
do
echo "doing salmon stuff..."
salmon quant -p 16 -i GRCg7b_transcripts.salmon-index-plus-decoys.idx -l ISF -1 $line\_R1.fastq.gz -2 $line\_R2.fastq.gz --softclip --writeMappings -g GCF_016699485.2_bGalGal1.mat.broiler.GRCg7b_genomic.gff3 -o $line\.quant | samtools sort -@16 -O BAM -o $line\.bam
done
