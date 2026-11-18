#!/usr/bin/env bash

#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=24
#SBATCH --mem=64GB
#SBATCH --partition=whelan,batch
#SBATCH --job-name=AL-pearlshell-nrule
#SBATCH --mail-type=ALL
#SBATCH --mail-user=Whelan.105@osu.edu
#SBATCH --hint=multithread

# Paris et al. (2017) n rule at the chosen M = 2: compare n = M - 1, M, M + 1.
# n = 2 is already 2_mmar / 2_mmar_R80, so only n = 1 and n = 3 are run here.
# -X "populations:-R 0.8" makes denovo_map.pl's own populations step Paris's r80 filter
# (loci in >= 80% of ALL samples), so 2_mmar_n1/ and 2_mmar_n3/ can be passed straight to
# paris_optimize.py. Do NOT add -r 0.8: populations would then apply both filters.
# Afterwards:
#   python paris_optimize.py {1..6}_mmar_R80 2_mmar_n1 2_mmar_n3 -o paris_results

# cd to dir from which I submit job
cd $SLURM_SUBMIT_DIR

for n in 1 3
do
mkdir 2_mmar_n${n}
denovo_map.pl --samples ./data --popmap ./popmap_mmar.txt -m 3 -M 2 -n ${n} -T 24 -o ./2_mmar_n${n} --paired -X "populations:-R 0.8"
done
