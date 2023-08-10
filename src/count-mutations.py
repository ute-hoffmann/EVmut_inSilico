#!/usr/bin/env python3
"""

"""

import os
import pandas as pd

path_file_1 = "results/Sequences_Fitness_Point-Mutations.csv"
data_WT_aa = {}
with open(path_file_1) as f:
	linecounter = 1
	for line in f:
		if linecounter == 1:
			linecounter += 1
			continue
		# structure of file: Fitness,NumberMutations,PointMutations,Sequence\n
		line_split = line.strip("\n").split(",")
		pointMutations = line_split[2]
		for mutation in pointMutations.split(" "):
			first_amino_acids = mutation.strip(" ")[:-1]
			mutated_to = mutation.strip(" ")[-1:]
			if first_amino_acids in data_WT_aa.keys():
				data_WT_aa[first_amino_acids][0] += 1
				if mutated_to not in data_WT_aa[first_amino_acids][1]:
					data_WT_aa[first_amino_acids][1].append(mutated_to)
			else:
				data_WT_aa[first_amino_acids] = [1, [mutated_to]]
		
data_WT_file = {"Mutation":[], "Number":[], "Mutated_to":[]}
for mut in data_WT_aa.keys():
	data_WT_file["Mutation"].append(mut)
	data_WT_file["Number"].append(str(data_WT_aa[mut][0]))
	string_mutations = ""
	for i in data_WT_aa[mut][1]:
		string_mutations = string_mutations + i + " "
	data_WT_file["Mutated_to"].append(string_mutations)

df = pd.DataFrame(data_WT_file)
df.to_csv('results/countMutations_limitedPos.csv', index=False)
