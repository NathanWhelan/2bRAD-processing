#!/usr/bin/env bash

#SBATCH --time=150:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=24
#SBATCH --mem=32GB
#SBATCH --partition=whelan,batch
#SBATCH --job-name=pars


# cd to dir from which I submit job
cd $SLURM_SUBMIT_DIR


# Assumes demultiplexed reads are in a folder called data.
# assumes popmap is called "popmap_params-test.txt"
for i in {1..6}
do
mkdir $i\_params-test
denovo_map.pl --samples ./data --popmap ./popmap_params-test.txt -m 3 -M $i -n $i -T 24 -o ./$i\_params-test -X populations:"-R 0.80"
done

