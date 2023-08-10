#!/usr/bin/env python3
"""
Read in data from replicates and store important information (replicate number, fitness values of accepted sequences, temperature) in dataframe, safe as csv

Author: Ute Hoffmann
"""

import os
import pandas as pd

path_directory = "results/replication_3mutations"
files = os.listdir(path_directory)

data_holder = {"Fitness":[], "Step":[], "Replicate":[], "Temperature":[]}

for file in files:
	fileName_split = file.split("_")
	replicate_number = fileName_split[0]
	temperature = fileName_split[1]
	with open(path_directory + "/" + file) as f:
		linecounter = 1
		for line in f:
			if linecounter == 1:
				linecounter += 1
				continue
			# structure of file: sequences\tScores\taccept\tstep
			line_split = line.strip("\n").split("\t")
			if line_split[2] == "True":
				score = line_split[1]
				step = line_split[-1]
				data_holder["Fitness"].append(score)
				data_holder["Step"].append(step)
				data_holder["Replicate"].append(replicate_number)
				data_holder["Temperature"].append(temperature)
				
df = pd.DataFrame(data_holder)
df.to_csv('results/2000replicates_500steps_3trust.csv', index=False)
