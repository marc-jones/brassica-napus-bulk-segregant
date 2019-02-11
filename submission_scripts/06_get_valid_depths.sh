#!/bin/bash

DARMOR='ETNV208_DSW57010_HFL2JCCXY_L7'
CABRIOLET='ETNV203_DSW57009_HFL2JCCXY_L6'

PILEUPDIR='./mpileup_unique_2018_04_20_15_02'

MPILEUP_REF=$(find $PILEUPDIR -name "mpileup_reference")

DARMOR_NUM=$(expr $(grep -n $DARMOR $MPILEUP_REF | cut -c1) - 1)
CABRIOLET_NUM=$(expr $(grep -n $CABRIOLET $MPILEUP_REF | cut -c1) - 1)

OUTPUT_DIR=$(date +depth_unique_%Y_%m_%d_%H_%M) ### Changed
mkdir -p ./$OUTPUT_DIR/run_logs
mkdir -p ./$OUTPUT_DIR/slurm_scripts
mkdir -p ./$OUTPUT_DIR/output_data

echo \
"#!/bin/bash -e
#SBATCH -p nbi-long
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem 32000
#SBATCH -o ./$OUTPUT_DIR/run_logs/nv5.%N.%j.out
#SBATCH -e ./$OUTPUT_DIR/run_logs/nv5.%N.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=marc.jones@jic.ac.uk

python pileup_valid_depth.py \
    ./$OUTPUT_DIR/output_data/nv5.csv \
    $DARMOR_NUM \
    $CABRIOLET_NUM \
    $(expr $(grep -n "ETNVLP6_DSW57008_HFL2JCCXY_L6" $MPILEUP_REF | cut -c1) - 1) \
    $(expr $(grep -n "ETNVEP5_DSW57007_HFL2JCCXY_L7" $MPILEUP_REF | cut -c1) - 1)

" > ./$OUTPUT_DIR/slurm_scripts/nv5.sh
JOBID=$(sbatch ./$OUTPUT_DIR/slurm_scripts/nv5.sh | \
    awk '{print($4)}')
echo "Submitted job $JOBID"

echo \
"#!/bin/bash -e
#SBATCH -p nbi-long
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem 32000
#SBATCH -o ./$OUTPUT_DIR/run_logs/v510.%N.%j.out
#SBATCH -e ./$OUTPUT_DIR/run_logs/v510.%N.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=marc.jones@jic.ac.uk

python pileup_valid_depth.py \
    ./$OUTPUT_DIR/output_data/v510.csv \
    $DARMOR_NUM \
    $CABRIOLET_NUM \
    $(expr $(grep -n "TD170810375_DSW50199_H7JCVCCXY_L3" $MPILEUP_REF | cut -c1) - 1) \
    $(expr $(grep -n "TD170810374_DSW50198_H7JCVCCXY_L2" $MPILEUP_REF | cut -c1) - 1)

" > ./$OUTPUT_DIR/slurm_scripts/v510.sh
JOBID=$(sbatch ./$OUTPUT_DIR/slurm_scripts/v510.sh | \
    awk '{print($4)}')
echo "Submitted job $JOBID"

echo \
"#!/bin/bash -e
#SBATCH -p nbi-long
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem 32000
#SBATCH -o ./$OUTPUT_DIR/run_logs/v5.%N.%j.out
#SBATCH -e ./$OUTPUT_DIR/run_logs/v5.%N.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=marc.jones@jic.ac.uk

python pileup_valid_depth.py \
    ./$OUTPUT_DIR/output_data/v5.csv \
    $DARMOR_NUM \
    $CABRIOLET_NUM \
    $(expr $(grep -n "TD170810376_DSW50200_H7JCVCCXY_L3" $MPILEUP_REF | cut -c1) - 1) \
    $(expr $(grep -n "TD170810373_DSW50197_H7JCVCCXY_L2" $MPILEUP_REF | cut -c1) - 1)

" > ./$OUTPUT_DIR/slurm_scripts/v5.sh
JOBID=$(sbatch ./$OUTPUT_DIR/slurm_scripts/v5.sh | \
    awk '{print($4)}')
echo "Submitted job $JOBID"
