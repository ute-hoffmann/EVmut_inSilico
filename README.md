# EVmutation-based in silico evolution

Code in this repository helps to identify a zero-shot prediction ([Wittmann et al., 2021](https://doi.org/10.1016/j.cels.2021.07.008)) for a higher-order mutant variant library in case there is no prior experimental information. To do so, it uses a Metropolis-Hastings Markov algorithm to navigate a fitness landscape given by the phylogeny-guided model calculated by [EVcouplings](https://v2.evcouplings.org/). 

# Installation

Requirements: bash, python3/miniconda, R

Clone this repository

```bash
git clone https://github.com/ute-hoffmann/EVmut_inSilico
```

Set up conda environment, run this command on shell from within folder in which file evmut_insilico.yml is located

```bash
conda env create -f evmut_insilico.yml
```

# Use

For a new protein sequence: retrieve amino acid sequence, run EVcouplings web server https://v2.evcouplings.org/
Take into account that enough sequences for multiple sequence alignment, probably not feasible for proteins with to few homologs

download model into input folder, adjust file iterate_inSilicoEvolution.sh to number of iterations, temperatures, output directory, where to find .model file, original amino acid sequence

Activate conda environment 

```bash
conda activate evmut_insilico
bash iterate_insilicoEvolution.sh
```

Plotting of fitness values for different fitness values: 
prepare_df_for_plotting.py

then plot_fitness_perStep.R
run plot_fitness_perStep.R

choose temperature for MCMC sampling (not too much noise, but big enough steps possible)

use extract_best_in_trajectory.py to extract best sequences in trajectory

use extract_point-mutations_number-mutations.py to extract point mutations of these

use count_mutations.py to count how often which mutation showed up


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
