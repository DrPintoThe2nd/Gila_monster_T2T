#!/bin/bash

#SBATCH -c 24                  		# number of cores
#SBATCH --mem=48GB              	# specify additional memory
#SBATCH -t 1-00:00:00          	    # wall time (D-HH:MM)
#SBATCH -p public      	        # Uses the
#SBATCH -q public     	            # Submits job
#SBATCH -o minimap.%j.out       	 	# STDOUT (%j = JobId)
#SBATCH -e minimap.%j.err       	 	# STDERR (%j = JobId)
#SBATCH --mail-type=ALL             # Send a notification when a job starts, stops, or fails
#SBATCH --mail-user=cpwebste@asu.edu # send-to address
#SBATCH --export=NONE

# define variables
input_ref='/scratch/bpinto2/Carmen_share/Heloderma_suspectum.final.reference.fasta'

input_hifi_reads='/scratch/bpinto2/Carmen_share/Gila_monster_HiFi.fastq.gz'
output_hifi_bam='/scratch/cpwebste/gila_monster/gila_final_ref_HiFi.bam'

input_ont_reads='/scratch/bpinto2/Carmen_share/Gila_monster.combined.UL.5000.simplex.fq.gz'
output_ont_bam='/scratch/cpwebste/gila_monster/gila_final_ref_ONT.bam'

module purge
module load mamba/latest

source activate nanopolish
# conda install bioconda::sambamba
# conda install bioconda::minimap2=2.30
# conda install bioconda::samtools=1.22.1

# PB-HiFi alignment
minimap2 -t 24 -ax map-hifi --eqx --secondary=no ${input_ref} ${input_hifi_reads} |\
       	samtools sort -@ 24 -O BAM -o ${output_hifi_bam}

# ONT-UL alignment
minimap2 -t 24 -ax lr:hq --eqx --secondary=no ${input_ref} ${input_ont_reads} |\
       	samtools sort -@ 24 -O BAM -o ${output_ont_bam}
# Nanopore Q20 genomic reads (v2.27 or later)

# index PB-HiFi alignment
sambamba index -p -t4 ${output_hifi_bam}
# index ONT-UL alignment
sambamba index -p -t4 ${output_ont_bam}

source deactivate

