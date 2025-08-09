#!/bin/bash

#SBATCH --account=bgmp
#SBATCH --partition=bgmp
#SBATCH --cpus-per-task=1
#SBATCH --mem=32GB
#SBATCH --mail-user=vini@uoregon.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name=demultiplex
#SBATCH --output=demultiplex_%j.out
#SBATCH --error=demultiplex_%j.err

/usr/bin/time -v \
./demultiplex.py \
-r1 '/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R1_001.fastq.gz' \
-r2 '/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R2_001.fastq.gz' \
-r3 '/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R3_001.fastq.gz' \
-r4 '/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R3_001.fastq.gz' \
-b '/projects/bgmp/shared/2017_sequencing/indexes.txt' \
-o 'real_run'