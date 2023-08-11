# EVmutation-based in silico evolution

Code in this repository helps to identify a zero-shot prediction ([Wittmann et al., 2021](https://doi.org/10.1016/j.cels.2021.07.008)) for a higher-order mutant variant library in case there is no prior experimental information. To do so, it uses a Metropolis-Hastings Markov algorithm to navigate a fitness landscape given by the phylogeny-guided model calculated by [EVcouplings](https://v2.evcouplings.org/). 

# System requirements

## Hardware

ran on normal laptop

## Software

Requirements: bash, Miniconda (or Anaconda), preferentially mamba, for some code R, R library tidyverse, required python libraries will be installed using conda environment

| Software | Version | Tested version | Note |
| -------- | ------- | -------------- | --------- |
| OS | | Ubuntu 22.04? | |
| Bash | 4.0 | ? | |
| Python | >3.6 | 3.6.2?, 3.9.2 ? | |
| R | 3.6.3 | 3.6.3, 4.1.1 | |
| GNU parallel | 20161222 | 20161222 | |


# Installation

Clone this repository

```bash
git clone https://github.com/ute-hoffmann/EVmut_inSilico
```

Set up conda environment, run this command on shell from within folder in which file evmut_insilico.yml is located - preferentially mamba to resolve environment!

```bash
mamba env create -f evmut_insilico.yml
```

# Use

For a new protein sequence: retrieve amino acid sequence, run EVcouplings web server https://v2.evcouplings.org/
Take into account that enough sequences for multiple sequence alignment, probably not feasible for proteins with to few homologs

download model into input folder, adjust file iterate_inSilicoEvolution.sh to number of iterations, temperatures, output directory, where to find .model file, original amino acid sequence

Activate conda environment, run from base folder

```bash
conda activate evmut_insilico
bash iterate_insilicoEvolution.sh
```

Plotting of fitness values for different fitness values: 
-j gives job_specific_name (in directory path!)
```bash
python src/prepare_df_for_plotting.py -j Gallionella_cbbM
```
creates file collect_data_plotting.csv with all data for replicates, steps, temperatures etc. 

required package: tidyverse
script assumes 500 steps (line 13): adjust accordingly! takes a while!
```bash
Rscript src/plot_fitness_perStep.R Gallionella_cbbM
```

choose temperature for MCMC sampling (not too much noise, but big enough steps possible), here: T=0.03

![alt text](results/Gallionella_cbbM/collect_data_plotting.png "Average fitness of variants after certain number of steps in MCMC, for different temperatures")

--> implement: status updates! + better output names

```bash
python src/extract_and_count.py -j Gallionella_cbbM -s input/Gallionella_Rubisco.txt -T 0.03
```

Run in silico evolution again with only limited number of amino acids, check which sequences pop up most frequently together and check their fitness score, A230E does at least pop up among most highly scoring sequences

## Output

Files: 
Trajectories
file with best in trajectories
file with amino acid substitutions present in best in trajectories
file with how frequent these substitutions showed up
png file with step ~ fitness

# Citing supporting repositories

EVcouplings (evcouplings.couplings is used) and also jax_unirep - Unirep repository

Loosely based on [github.com/Asplund-Samuelsson/furee/blob/master/source/in_silico_evolution.py](https://github.com/Asplund-Samuelsson/furee/blob/master/source/in_silico_evolution.py)
EVmutation code from [github.com/fhalab/MLDE/blob/39e9edccc346119a834c62677d16e39dd49dfbdb/code/zero_shot/zero_shot_predictor.py#L29](https://github.com/fhalab/MLDE/blob/39e9edccc346119a834c62677d16e39dd49dfbdb/code/zero_shot/zero_shot_predictor.py#L29)

MCMC code taken from [github.com/ElArkk/jax-unirep/tree/master](https://github.com/ElArkk/jax-unirep/tree/master)
