#!/bin/bash

REFERENCE='/usr/users/cbu/jonesd/tophat/brassica_napus/genome/b_napus'

### Create an output directory using the current date and time
OUTPUT_DIR=$(date +mpileup_unique_%Y_%m_%d_%H_%M)
mkdir -p ./$OUTPUT_DIR/run_logs
mkdir -p ./$OUTPUT_DIR/slurm_scripts
mkdir -p ./$OUTPUT_DIR/output_data

BAM_FILES=$(find ./ -name "*_unique.bam")

find ./ -name "*_unique.bam" > $OUTPUT_DIR/mpileup_reference

echo \
"#!/bin/bash -e
#SBATCH -p nbi-long
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 1
#SBATCH --mem 32000
#SBATCH -o ./$OUTPUT_DIR/run_logs/mpileup_unique.%N.%j.out
#SBATCH -e ./$OUTPUT_DIR/run_logs/mpileup_unique.%N.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=marc.jones@jic.ac.uk

source samtools-1.5

# perform a scan for varients
srun samtools mpileup -f $REFERENCE.fa \
    $(echo $BAM_FILES) > ./$OUTPUT_DIR/output_data/combined.pileup

" > ./$OUTPUT_DIR/slurm_scripts/mpileup_unique.sh
JOBID=$(sbatch ./$OUTPUT_DIR/slurm_scripts/mpileup_unique.sh | \
    awk '{print($4)}')
echo "Submitted job $JOBID"
