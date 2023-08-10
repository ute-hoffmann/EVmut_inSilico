#!/bin/bash

# would be nice to have: read in arguments from command line

# Create output directory
mkdir -p results/replication_3mutations

# Iterate over 2001 replicates
for i in {0..2000}; do
  # Iterate over three temperatures
  for T in 0.3 0.03 0.003; do
      # Evolve based on difference of candidate and best sequence scores
      Outfile="results/replication_3mutations/${i}_${T}_Difference.tab"
      python in_silico_evolution_EVcouplings.py \
        -s 500 -t 2 -T $T \ # -s : steps, -t : trust radius, -T : temperature
        -EV ../input/EVmut_output/92cbb84d7f8a4239b7fb4f195c8bc187_TARGET_b0.5.model \ # model calculated by web server https://v2.evcouplings.org/
        Gallionella_Rubisco.txt \ # original amino acid sequence
        $Outfile
  done
done

