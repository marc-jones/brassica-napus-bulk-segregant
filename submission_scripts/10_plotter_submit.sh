#!/bin/bash -e
#SBATCH -p nbi-medium
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 1
#SBATCH --mem 32000
#SBATCH -o ./plotting/plotting.%N.%j.out
#SBATCH -e ./plotting/plotting.%N.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=marc.jones@jic.ac.uk

source R-3.4.0-20170516

Rscript plotter.R
