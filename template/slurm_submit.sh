#!/bin/bash
#SBATCH --job-name=RNA_seq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=40g
#SBATCH --time=12:00:00
#SBATCH --output=RNA_seq.out
#SBATCH --error=RNA_seq.err

WORKDIR=$SLURM_SUBMIT_DIR

# Arguments
PIPELINE=$(grep -Po '"pipeline": *\K"[^"]*' snake_conf.json | tr -d '"')
RUN_TYPE=$(grep -Po '"run_type": *\K"[^"]*' snake_conf.json | tr -d '"')

## Set working directory to scratch
SCRATCHDIR=/mnt/hpcscratch/cgonzalez/R24/DICE_GALAXY/remap/BEN_RNAseq

# Create scratch directory if it does not exist
mkdir -p $SCRATCHDIR

# Copy the pipeline to scratch
cp -r $PIPELINE/scripts $SCRATCHDIR/
cp $PIPELINE/Snakefile_$RUN_TYPE $SCRATCHDIR/
cp snake_conf.json $SCRATCHDIR/
cp chr_name_conv.txt $SCRATCHDIR/
ln -s $WORKDIR/1.Fastq_input $SCRATCHDIR/

# Change to scratch directory
cd $SCRATCHDIR

# Remove necessary to make sure the run info is correct
rm -rf logs/*
mkdir -p logs

python scripts/move_fastq.py -json snake_conf.json

snakemake --profile /mnt/BioAdHoc/Groups/vd-vijay/Cristian/DICE_GALAXY/RNA_seq/profiles/slurm --configfile snake_conf.json --snakefile $PIPELINE/Snakefile_$RUN_TYPE --stats logs/snakemake.stats >& logs/snakemake.log 

rm -rf scripts
