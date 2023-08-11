#!/bin/bash

# adjust paths to files and settings here
EVmodel="EVmut_output/92cbb84d7f8a4239b7fb4f195c8bc187_TARGET_b0.5.model" # model calculated by web server https://v2.evcouplings.org/, is assumed to be located in input directory
input_sequence="Gallionella_Rubisco.txt" # original amino acid sequence, not with fasta header or similar, is assumed to be located in input directory
number_replicates=2000 # how many trajectories should be run?
steps=500 # number of steps of MCMC chain
trust=2 # how many amino acid exchanges may be introduced into sequence
job_specific_name="Gallionella_cbbM" # name included in folder structure

# Create output directory
mkdir -p results/${job_specific_name}/trajectories

# Iterate over certain number of replicates
for((i=1; i<=$number_replicates; i++)); do
  # Iterate over different temperatures (include or exclude from for loop if needed)
  for T in 0.3 0.03 0.003; do
      Outfile="results/"$job_specific_name"/trajectories/${i}_${T}_Difference.tab"
      python src/in_silico_evolution_EVcouplings.py -s $steps -t $trust -T $T -EV input/$EVmodel input/$input_sequence $Outfile
  done
done

