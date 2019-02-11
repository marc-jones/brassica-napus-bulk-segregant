#!/bin/bash

TOPHAT_DIR='/usr/users/cbu/jonesd/tophat'
IRWIN_SHARE='/nbi/group-data/ifs/JIC/Judith-Irwin'
REFERENCE='/usr/users/cbu/jonesd/tophat/brassica_napus/genome/b_napus'
READ_DIR_1=$IRWIN_SHARE'/2017_11_21_eleri_sequencing/C101HW17071262/raw_data'
READ_DIR_2=$IRWIN_SHARE'/Judith/C101HW17071262/raw_data/'
VARSCAN_DIR=$IRWIN_SHARE'/2017_11_21_eleri_sequencing'
CORES=32

### Create an output directory using the current date and time
OUTPUT_DIR=$(date +out_%Y_%m_%d_%H_%M)
mkdir -p ./$OUTPUT_DIR/run_logs
mkdir -p ./$OUTPUT_DIR/slurm_scripts
mkdir -p ./$OUTPUT_DIR/output_data

UNIQUE_IDENTIFIERS=($(find $READ_DIR_1 $READ_DIR_2 -name "*.fq.gz" | \
                      awk '{split($0, a, "/");
                          gsub(/_[1-2]\.fq\.gz/, "", a[length(a)]);
                          print(a[length(a)])}' | \
                      sort -u))

for (( ident_idx=0; ident_idx<${#UNIQUE_IDENTIFIERS[@]}; ident_idx++ ));
do
    IDENTIFIER=${UNIQUE_IDENTIFIERS[$ident_idx]}

    READS_1=`find $READ_DIR_1 $READ_DIR_2 -regex ".*"$IDENTIFIER"_1\.fq\.gz"`
    READS_2=`find $READ_DIR_1 $READ_DIR_2 -regex ".*"$IDENTIFIER"_2\.fq\.gz"`

    WORKING_DIR='./'$OUTPUT_DIR'/output_data/'$IDENTIFIER

    mkdir -p $WORKING_DIR

    echo \
"#!/bin/bash -e
#SBATCH -p nbi-medium
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c $CORES
#SBATCH --mem 32000
#SBATCH -o ./$OUTPUT_DIR/run_logs/$IDENTIFIER.%N.%j.out
#SBATCH -e ./$OUTPUT_DIR/run_logs/$IDENTIFIER.%N.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=marc.jones@jic.ac.uk

export PATH=$TOPHAT_DIR/bin/bowtie2-2.2.3:$PATH
source hpccore-5
source samtools-1.5
source jdk-1.7.0_25

cd $WORKING_DIR

# run bowtie
srun bowtie2 -p $CORES \
    --rg-id $IDENTIFIER \
    --rg SM:$IDENTIFIER \
    -x $REFERENCE \
    -1 $READS_1 \
    -2 $READS_2 \
    -S $IDENTIFIER.sam

# convert from sam to bam
srun samtools view -b -@ $CORES $IDENTIFIER.sam > $IDENTIFIER.bam

# sort the bam file
srun samtools sort -o $IDENTIFIER.sorted.bam -@ $CORES $IDENTIFIER.bam

# delete the unnecessary sam and bam files
rm -f $IDENTIFIER.sam $IDENTIFIER.bam

# index the sorted bam file
srun samtools index $IDENTIFIER.sorted.bam $IDENTIFIER.sorted.bai

# calculate some stats on the alignment
srun samtools flagstat $IDENTIFIER.sorted.bam > $IDENTIFIER\_stats.out
srun samtools depth $IDENTIFIER.sorted.bam > $IDENTIFIER\_depth.out

# perform a scan for varients
srun samtools mpileup -f $REFERENCE.fa \
    $IDENTIFIER.sorted.bam > $IDENTIFIER.pileup

    " > ./$OUTPUT_DIR/slurm_scripts/$IDENTIFIER.sh
    JOBID=$(sbatch ./$OUTPUT_DIR/slurm_scripts/$IDENTIFIER.sh | \
        awk '{print($4)}')
    echo "Submitted job $JOBID"

done
