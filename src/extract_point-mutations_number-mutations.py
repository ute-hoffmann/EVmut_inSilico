#!/usr/bin/env python3
"""
Read in file with best trajectories and return file with point mutations, number of mutations, fitness score, sequence
"""

import os
import pandas as pd

path_file = "results/Temperature_0-03_bestSequencesTrajectory.csv"

# Load target sequence
sequence_file = "Gallionella_Rubisco.txt"
with open(sequence_file) as s:
    starting_sequence = s.read().strip()

data_holder = {"Fitness":[], "NumberMutations":[], "PointMutations":[], "Sequence":[]}

with open(path_file) as f:
	linecounter = 1
	for line in f:
		if linecounter == 1:
			linecounter += 1
			continue
		# structure of file: Fitness,Step,Replicate,Sequence\n
		line_split = line.strip("\n").split(",")
		score = float(line_split[0])
		sequence=line_split[-1]
		mutation_list = ""
		mutation_counter = 0
		for i in range(len(starting_sequence)):
        		if starting_sequence[i] != sequence[i]:
              			mutation_list += " " + starting_sequence[i] + str(i+1) + sequence[i]
              			mutation_counter += 1
		data_holder["Fitness"].append(score)
		data_holder["Sequence"].append(sequence)
		data_holder["NumberMutations"].append(mutation_counter)
		data_holder["PointMutations"].append(mutation_list)	
	
df = pd.DataFrame(data_holder)
df.to_csv('results/2022-05-04_Sequences_Fitness_Point-Mutations.csv', index=False)
