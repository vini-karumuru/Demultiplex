#!/bin/bash
#SBATCH --account=bgmp
#SBATCH --partition=bgmp
#SBATCH --cpus-per-task=8
#SBATCH --mem=64GB
#SBATCH --mail-user=vini@uoregon.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name=plot_qs_dist
#SBATCH --output=plot_qs_dist_%j.out
#SBATCH --error=plot_qs_dist_%j.err

/usr/bin/time -v \
./plot_quality_distribution.py \
-i /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R1_001.fastq.gz \
-o R1_qs_dist.png

/usr/bin/time -v \
./plot_quality_distribution.py \
-i /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R2_001.fastq.gz \
-o R2_qs_dist.png

/usr/bin/time -v \
./plot_quality_distribution.py \
-i /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R3_001.fastq.gz \
-o R3_qs_dist.png

/usr/bin/time -v \
./plot_quality_distribution.py \
-i /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R4_001.fastq.gz \
-o R4_qs_dist.png