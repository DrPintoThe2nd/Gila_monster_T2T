#!/bin/bash

#SBATCH -c 12                  		# number of cores
#SBATCH --mem=12GB              	# specify additional memory to accomodate bbmap
#SBATCH -t 0-00:30:00                   # wall time (D-HH:MM)
#SBATCH -p public      	        # Uses the
#SBATCH -q public     	                # Submits job
#SBATCH -o tidk.%j.out      		# STDOUT (%j = JobId)
#SBATCH -e tidk.%j.err        	# STDERR (%j = JobId)
#SBATCH --mail-type=ALL                 # Send a notification when a job starts, stops, or fails
#SBATCH --mail-user=bjp004@morningside.edu
#SBATCH --export=NONE

module purge

source activate liftoff

tidk search --string CCCTAA --output Heloderma_suspectum.06-2026.hap1 --dir . Heloderma_suspectum.06-2026.hap1.fasta
tidk plot --tsv Heloderma_suspectum.06-2026.hap1_telomeric_repeat_windows.tsv --output Heloderma_suspectum.06-2026.hap1 --height 120 --width 800

tidk search --string CCCTAA --output Heloderma_suspectum.10-2025.hap2 --dir . Heloderma_suspectum_10-3-25.hap2.fasta
tidk plot --tsv Heloderma_suspectum.10-2025.hap2_telomeric_repeat_windows.tsv --output Heloderma_suspectum.10-2025.hap2 --height 120 --width 800
