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

#########do some quick QC
#fastqc -t 16 *.fastq.gz
#multiqc .

#########make informative transcriptomes for each species
#gffread -g Eublepharis_macularius_reference.genome.fasta -x Eublepharis_macularius_reference.gffread.fasta --no-pseudo Eublepharis_macularius_reference.genome.gff
#gffread -g Pantherophis_guttatus_reference.genome.fasta -x Pantherophis_guttatus_reference.gffread.fasta --no-pseudo Pantherophis_guttatus_reference.genome.gff
#gffread -g Python_molurus_reference.genome.fasta -x Python_molurus_reference.gffread.fasta --no-pseudo Python_molurus_reference.genome.gff
#gffread -g Vipera_berus_reference.genome.fasta -x Vipera_berus_reference.gffread.fasta --no-pseudo Vipera_berus_reference.genome.gff

#########make index for each species
#salmon index -p 16 -k 31 -t Eublepharis_macularius_reference.gffread.fasta -i Eublepharis_macularius_reference.gffread.idx
#salmon index -p 16 -k 31 -t Pantherophis_guttatus_reference.gffread.fasta -i Pantherophis_guttatus_reference.gffread.idx
#salmon index -p 16 -k 31 -t Python_molurus_reference.gffread.fasta -i Python_molurus_reference.gffread.idx
#salmon index -p 16 -k 31 -t Vipera_berus_reference.gffread.fasta -i Vipera_berus_reference.gffread.idx

###############calc salmon Echis_coloratus using Vipera_berus_reference
ls ERR*_Echis_coloratus_*_R1.fastq.gz > Echis_coloratus.txt
sed -i 's/_R1.fastq.gz//g' Echis_coloratus.txt
sed -i 's/@//' Echis_coloratus.txt

cat Echis_coloratus.txt | while read line
do
echo "doing salmon stuff..."
salmon quant -p 16 -i Vipera_berus_reference.gffread.idx -l A -1 $line\_R1.fastq.gz -2 $line\_R2.fastq.gz --softclip --writeMappings -g Vipera_berus_reference.genome.gff -o $line\.quant | samtools sort -@16 -O BAM -o $line\.bam
done

###############calc salmon Eublepharis_macularius using Eublepharis_macularius_reference
ls ERR*_Eublepharis_macularius_*_R1.fastq.gz > Eublepharis_macularius.txt
sed -i 's/_R1.fastq.gz//g' Eublepharis_macularius.txt
sed -i 's/@//' Eublepharis_macularius.txt

cat Eublepharis_macularius.txt | while read line
do
echo "doing salmon stuff..."
salmon quant -p 16 -i Eublepharis_macularius_reference.gffread.idx -l A -1 $line\_R1.fastq.gz -2 $line\_R2.fastq.gz --softclip --writeMappings -g Eublepharis_macularius_reference.genome.gff -o $line\.quant | samtools sort -@16 -O BAM -o $line\.bam
done

###############calc salmon Opheodrys_aestivus using Pantherophis_guttatus_reference
ls ERR*_Opheodrys_aestivus_*_R1.fastq.gz > Opheodrys_aestivus.txt
sed -i 's/_R1.fastq.gz//g' Opheodrys_aestivus.txt
sed -i 's/@//' Opheodrys_aestivus.txt

cat Opheodrys_aestivus.txt | while read line
do
echo "doing salmon stuff..."
salmon quant -p 16 -i Pantherophis_guttatus_reference.gffread.idx -l A -1 $line\_R1.fastq.gz -2 $line\_R2.fastq.gz --softclip --writeMappings -g Pantherophis_guttatus_reference.genome.gff -o $line\.quant | samtools sort -@16 -O BAM -o $line\.bam
done

###############calc salmon Pantherophis_guttatus using Pantherophis_guttatus_reference
ls ERR*_Pantherophis_guttatus_*_R1.fastq.gz > Pantherophis_guttatus.txt
sed -i 's/_R1.fastq.gz//g' Pantherophis_guttatus.txt
sed -i 's/@//' Pantherophis_guttatus.txt

cat Pantherophis_guttatus.txt | while read line
do
echo "doing salmon stuff..."
salmon quant -p 16 -i Pantherophis_guttatus_reference.gffread.idx -l A -1 $line\_R1.fastq.gz -2 $line\_R2.fastq.gz --softclip --writeMappings -g Pantherophis_guttatus_reference.genome.gff -o $line\.quant | samtools sort -@16 -O BAM -o $line\.bam
done

###############calc salmon Python_regius using Python_molurus_reference
ls ERR*_Python_regius_*_R1.fastq.gz > Python_regius.txt
sed -i 's/_R1.fastq.gz//g' Python_regius.txt
sed -i 's/@//' Python_regius.txt

cat Python_regius.txt | while read line
do
echo "doing salmon stuff..."
salmon quant -p 16 -i Python_molurus_reference.gffread.idx -l A -1 $line\_R1.fastq.gz -2 $line\_R2.fastq.gz --softclip --writeMappings -g Python_molurus_reference.genome.gff -o $line\.quant | samtools sort -@16 -O BAM -o $line\.bam
done

