#!/bin/bash

# adjust paths to files and settings
EVmodel="EVmut_output/92cbb84d7f8a4239b7fb4f195c8bc187_TARGET_b0.5.model" # model calculated by web server https://v2.evcouplings.org/, is assumed to be located in input directory
input_sequence="Gallionella_Rubisco.txt" # original amino acid sequence, not with fasta header or similar, is assumed to be located in input directory
file_positions_interest="input/positions_of_interest.txt"

number_replicates=750 # how many trajectories should be run?
job_specific_name="Gallionella_cbbM_limitedPos"
temperature=0.03
steps=500
trust=3

# Create output directory
mkdir -p results/${job_specific_name}/trajectories

# Iterate over certain number of replicates
for((i=1; i<=$number_replicates; i++)); do
  # Iterate over one temperature (which is technically not iterating then, I guess)
  for T in $temperature; do
      # Evolve based on difference of candidate and best sequence scores
      Outfile="results/"$job_specific_name"/trajectories/${i}_${T}_Difference.tab"
      python src/in_silico_evolution_zeroShot_EVcouplings_limitedPositions.py \
        -s $steps -t $trust -T $T \
        -EV input/$EVmodel \
        -i input/$input_sequence \
        -o $Outfile -posInt $file_positions_interest
  done
done

