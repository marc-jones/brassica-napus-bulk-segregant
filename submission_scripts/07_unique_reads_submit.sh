#!/bin/bash

OUTPUT_DIR=$(date +unique_%Y_%m_%d_%H_%M)

mkdir -p ./$OUTPUT_DIR/run_logs
mkdir -p ./$OUTPUT_DIR/slurm_scripts

BAM_FILES=($(find ./ -name *.bam))

# for (( i=0; i<1; i++ ))
for (( i=0; i<${#BAM_FILES[@]}; i++ ))
do
    BASE=$(basename ${BAM_FILES[$i]})
    BAM_FILE=${BAM_FILES[$i]}
    echo \
"#!/bin/bash -e
#SBATCH -p nbi-long
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 20
#SBATCH --mem 32000
#SBATCH -o ./$OUTPUT_DIR/run_logs/${BASE:0:8}.%N.%j.out
#SBATCH -e ./$OUTPUT_DIR/run_logs/${BASE:0:8}.%N.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=marc.jones@jic.ac.uk

source samtools-1.6

srun samtools view -b -q 42 -@ 20 ${BAM_FILES[$i]} | \
    samtools sort -o ${BAM_FILE/\.bam/_unique.bam} -@ 20 -

    " > ./$OUTPUT_DIR/slurm_scripts/${BASE:0:8}_filter_bam.sh

    sbatch ./$OUTPUT_DIR/slurm_scripts/${BASE:0:8}_filter_bam.sh
done
