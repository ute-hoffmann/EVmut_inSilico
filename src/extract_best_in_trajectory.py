#!/usr/bin/env python3
"""
Read in data for T=0.03, extract best protein seqeuence per trajectory, save data frame with sequence, fitness, replicate number, step
"""

import os
import pandas as pd

path_directory = "results/replication_3mutations"
files = os.listdir(path_directory)

data_holder = {"Fitness":[], "Step":[], "Replicate":[], "Sequence":[]}

for file in files:
	fileName_split = file.split("_")
	replicate_number = fileName_split[0]
	if fileName_split[1] != "0.03":
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
df.to_csv('results/Temperature_0-03_bestSequencesTrajectory.csv', index=False)
