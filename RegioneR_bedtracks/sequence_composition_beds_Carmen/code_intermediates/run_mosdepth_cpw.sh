#!/bin/bash

#SBATCH -c 18                       # number of cores
#SBATCH -N 1                        # number of Nodes
#SBATCH --mem=100G                 # specify additional memory to accomodate analysis
#SBATCH -t 0-02:00                  # wall time (D-HH:MM)
#SBATCH -p htc                      # Uses the high throughput partition 
#SBATCH -q public                   # Submits job
#SBATCH -o mosdepth.%j.out           # STDOUT (%j = JobId)
#SBATCH -e mosdepth.%j.err           # STDERR (%j = JobId)
#SBATCH --mail-type=ALL             # Send a notification when a job starts, stops, or fails
#SBATCH --mail-user=cpwebste@asu.edu # send-to address
#SBATCH --export=NONE
#SBATCH --job-name=mosdepth

### NAME VARIABLES
input_HiFi_bam='/scratch/cpwebste/gila_monster/gila_final_ref_HiFi.bam'
output_HiFi_pre='/scratch/cpwebste/gila_monster/mosdepth/PB_HiFi_gila_depth'

input_ONT_bam='/scratch/cpwebste/gila_monster/gila_final_ref_ONT.bam'
output_ONT_pre='/scratch/cpwebste/gila_monster/mosdepth/ONT_UL_gila_depth'

module purge
module load mamba/latest
echo "loaded mamba"
source activate mosdepth
echo "source activated"

### PB-HiFi
'''
# Use mosdepth to calculate seq depth at window size 100,000bp
MOSDEPTH_PRECISION=5 mosdepth -x -n \
	-t ${SLURM_CPUS_PER_TASK} \
	-Q 3 -b 100000 \
	${output_HiFi_pre}_100k \
       	${input_HiFi_bam}

# Use mosdepth to calculate seq depth at window size 500,000bp
MOSDEPTH_PRECISION=5 mosdepth -x -n \
	-t ${SLURM_CPUS_PER_TASK} \
	-Q 3 -b 500000 \
	${output_HiFi_pre}_500k \
       	${input_HiFi_bam}
'''
### ONT-UL

# Use mosdepth to calculate seq depth at window size 100,000bp
MOSDEPTH_PRECISION=5 mosdepth -x -n \
	-t ${SLURM_CPUS_PER_TASK} \
	-Q 3 -b 100000 \
	${output_ONT_pre}_100k \
       	${input_ONT_bam}

# Use mosdepth to calculate seq depth at window size 500,000bp
MOSDEPTH_PRECISION=5 mosdepth -x -n \
	-t ${SLURM_CPUS_PER_TASK} \
	-Q 3 -b 500000 \
	${output_ONT_pre}_500k \
       	${input_ONT_bam}

source deactivate
