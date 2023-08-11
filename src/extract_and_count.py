#!/usr/bin/env python3
"""
Read in data for T=0.03, extract best protein seqeuence per trajectory, save data frame with sequence, fitness, replicate number, step
"""

import os
import pandas as pd
import argparse

# Read arguments from the commandline
parser = argparse.ArgumentParser()

# Required input: Infile and project output (positional)
parser.add_argument(
    '-j', '--job_specific_name', type=str,
    help='Job specific name which is part of directory path.'
)

parser.add_argument(
    '-s', '--sequence_file', type=str,
    help='Text file with original amino acid sequence, no fasta header.'
)

parser.add_argument(
    '-T', '--temperature', type=str, default="0.03",
    help='MCMC temperature.'
)


# Parse arguments
args = parser.parse_args()

job_name = args.job_specific_name
sequence = args.sequence_file
temp = args.temperature

path_directory = "results/" + job_name + "/trajectories"
files = os.listdir(path_directory)

data_holder = {"Fitness":[], "Step":[], "Replicate":[], "Sequence":[]}

for file in files:
	fileName_split = file.split("_")
	replicate_number = fileName_split[0]
	if fileName_split[1] != temp:
		continue
	with open(path_directory + "/" + file) as f:
		linecounter = 1
		best_sequence = ""
		highscore = 0.0
		step = 0
		for line in f:
			if linecounter == 1:
				linecounter += 1
				continue
			# structure of file: sequences\tScores\taccept\tstep
			line_split = line.strip("\n").split("\t")
			score = float(line_split[1])
			if score < highscore:
				continue
			else:
				best_sequence = line_split[0]
				step = line_split[-1]
				highscore = score
	data_holder["Fitness"].append(highscore)
	data_holder["Step"].append(step)
	data_holder["Replicate"].append(replicate_number)
	data_holder["Sequence"].append(best_sequence)	
	
df = pd.DataFrame(data_holder)
df.to_csv("results/" + job_name + "/Temperature_"+temp +"_bestSequencesTrajectory.csv", index=False)

# Load target sequence
sequence_file = sequence
with open(sequence_file) as s:
    starting_sequence = s.read().strip()

# determine and count mutations in new sequences
data_holder_mut = {"Fitness":[], "NumberMutations":[], "PointMutations":[], "Sequence":[]}

for entry in range(len(data_holder["Fitness"])):
	score = data_holder["Fitness"][entry]
	sequence = data_holder["Sequence"][entry]
	mutation_list = ""
	mutation_counter = 0
	for i in range(len(starting_sequence)):
       		if starting_sequence[i] != sequence[i]:
       			mutation_list += " " + starting_sequence[i] + str(i+1) + sequence[i]
       			mutation_counter += 1
	data_holder_mut["Fitness"].append(score)
	data_holder_mut["Sequence"].append(sequence)
	data_holder_mut["NumberMutations"].append(mutation_counter)
	data_holder_mut["PointMutations"].append(mutation_list)	
	
df = pd.DataFrame(data_holder_mut)
df.to_csv("results/" + job_name + "/Sequences_Fitness_Point-Mutations.csv", index=False)

path_file_1 = "results/Gallionella_cbbM/Sequences_Fitness_Point-Mutations.csv"
data_WT_aa = {}

for entry in range(len(data_holder_mut["PointMutations"])):
	pointMutations = data_holder_mut["PointMutations"][entry]
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
df.to_csv("results/" + job_name + "/countMutations_limitedPos.csv", index=False)
