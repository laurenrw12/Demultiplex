#!/bin/bash

#SBATCH --account=bgmp
#SBATCH --partition=bgmp
#SBATCH -c 8

/usr/bin/time -v ./meandist.py \
-file /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R2_001.fastq.gz