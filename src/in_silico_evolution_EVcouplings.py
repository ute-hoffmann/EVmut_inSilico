#!/usr/bin/env python3
"""
Loosely based on https://github.com/Asplund-Samuelsson/furee/blob/master/source/in_silico_evolution.py
EVmutation code from https://github.com/fhalab/MLDE/blob/39e9edccc346119a834c62677d16e39dd49dfbdb/code/zero_shot/zero_shot_predictor.py#L29

Perform Metropolis-Hastings MCMC with EVcouplings model as scoring function to create zero-shot predictions (MCMC code taken from https://github.com/ElArkk/jax-unirep/tree/master)
i.e.: combine approach of Biswas et al., 2021 (Metropolis algorithm to perform in silico directed evolution) and Wittmann et al., 2021 (EVcouplings, zero-shot predictions)

author: Ute Hoffmann
"""

import numpy as np
import numpy.random as npr
import argparse
import pandas as pd
import os
from evcouplings.couplings import CouplingsModel
from collections import defaultdict
from typing import Callable, Dict
from tqdm.autonotebook import tqdm

# Read arguments from the commandline
parser = argparse.ArgumentParser()

# Required input: Infile and project output (positional)
parser.add_argument(
    'infile', type=str,
    help='Text file with one sequence on one line. Original fasta to start from.'
)

parser.add_argument(
    'outfile', type=str,
    help='Tab-delimited file with sampled sequences and scores.'
)

parser.add_argument(
    '-s', '--steps', type=int, default=10,
    help='Number of MCMC steps [10].'
)
parser.add_argument(
    '-t', '--trust', type=int, default=7,
    help='Trust radius [7].'
)
parser.add_argument(
    '-T', '--temperature', type=float, default=0.1,
    help='MCMC temperature; lower is slower [0.1].'
)
parser.add_argument(
    '-R', '--ratio', action='store_true', default=False,
    help='Use ratio instead of difference for sequence proposal rejection.'
)
parser.add_argument(
    '-EV', "--EVcouplings_model", type=str,
    help='Path to EVcouplings model (ends with .model), downloaded from EVcouplings webserver'
)

# Parse arguments
args = parser.parse_args()

sequence_file = args.infile
outfile = args.outfile
steps = args.steps
trust = args.trust
temperature = args.temperature
ratio = args.ratio
model_path = args.EVcouplings_model

# Load target sequence
with open(sequence_file) as s:
    starting_sequence = s.read().strip()

# Load a model for EVmutation
model = CouplingsModel(model_path)

# for EVmutation, only positions which were part of model may be proposed
letters_sorted = ""
for i in model.alphabet_map.keys():
    letters_sorted = letters_sorted + i
letters_sorted = sorted(letters_sorted)
aa_dict = model.alphabet_map

positions = model.index_list

def hamming_distance(s1: str, s2: str):
    """copied from jax_unirep, Return hamming distance between two strings of the same length."""
    return sum(c1 != c2 for c1, c2 in zip(s1, s2))

def propose(sequence: str): # patch in a way that sequences with only allowed indices and alphabet are used compared to jax_unirep version
    """
    Given a string, return a proposed mutated string.
    The proposed mutant is generated as follows:
    - Pick a position given the pick probabilities for each position
    - Given that position,
    pick a letter given the pick probabilities for each letter at that position
    that is not the existing letter.
    If no probabilites are specified for the positions or the pwm (position weight matrix),
    the sampler defaults to picking from an uniform distribution in both cases.
    :param sequence: The sequence to propose a new mutation on.
    :param pos_prob: Pick probability for each position in the sequence.
        Needs to be of shape (len(sequence), ),
        and probabilities need to sum to 1.
    :param pwm: Pick probability for each AA at each position.
        Needs to be of shape (len(sequence), 20),
        and probabilities at each position need to sum to 1.
        AA's need to be sorted alphabetically.
    :returns: A string.
    """
    if len(sequence) == 0:
        raise ValueError(
            "sequence passed into `propose` must be at least length 1."
        )

    # define uniform probability for positions and letters
    pos_prob = np.array([1.0/len(positions)] * len(positions))
    pwm = np.tile(np.array([[0.05] * 20]), (len(sequence), 1))

    position = positions[np.argmax(npr.multinomial(1, pos_prob))]
    new_sequence = ""
    for i, letter in enumerate(sequence):
        if (position-1) != i:
            new_sequence += letter
        else:
            letter_idx = aa_dict[letter]
            pwm[i, letter_idx] = 0
            new_letter_idx = np.argmax(npr.multinomial(1, pwm[i, :]))
            new_letter = letters_sorted[new_letter_idx]
            new_sequence += new_letter
    return new_sequence



    # Define sequence scoring function
def scoring_func(sequence: str):
    # create mutation list
    mutation_list = []
    for i in range(len(starting_sequence)):
        if starting_sequence[i] != sequence[i]:
              mutation_list.append((i+1, starting_sequence[i], sequence[i]))

    # Make the prediction, model loaded outside of function
    if len(mutation_list) > 0:
        delta_E, _, _ = model.delta_hamiltonian(mutation_list)
    else:
        delta_E = 0

    return delta_E

def is_accepted(best: float, candidate: float, temperature: float) -> bool:
    """
    from jax_unirep, uses differences and not ratios...
    Return boolean decision on whether the candidate mutant is accepted or not.

    This function checks whether we want to
    accept a new mutant proposed by our MMC sampler,
    by comparing its predicted activity
    to the current best mutants activity.

    :param best: Predicted activity of current best mutant
    :param candidate: Predicted activity of new candidate mutant
    :param temperature: Boltzmann distribution temperature.
        Controls acceptance probability.
        Low T decreases acceptance probability.
        High T increases acceptance probability.
    :returns bool: Whether or not candidate mutant was accepted
    """
    if (candidate-best)/temperature > 1:
        return True
    
    c = np.exp((candidate - best) / temperature)

    if c > 1:
        return True
    else:
        p = np.random.uniform(0, 1)
        if c >= p:
            return True
        else:
            return False

def sample_one_chain(
    starter_sequence: str,
    n_steps: int,
    scoring_func: Callable,
    is_accepted_kwargs: Dict = {},
    trust_radius: int = 7,
    propose_kwargs: Dict = {},
) -> Dict:
    """
    from repository: https://github.com/ElArkk/jax-unirep/tree/master
    Return one chain of MCMC samples of new sequences.

    Given a `starter_sequence`,
    this function will sample one chain of protein sequences,
    scored using a user-provided `scoring_func`.

    Design choices made here include the following.

    Firstly, we record all sequences that were sampled,
    and not just the accepted ones.
    This behaviour differs from other MCMC samplers
    that record only the accepted values.
    We do this just in case sequences that are still "good"
    (but not better than current) are rejected.
    The effect here is that we get a cluster of sequences
    that are one-apart from newly accepted sequences.

    Secondly, we check the Hamming distance
    between the newly proposed sequences and the original.
    This corresponds to the "trust radius"
    specified in the [jax-unirep paper](https://doi.org/10.1101/2020.01.23.917682).
    If the hamming distance > trust radius,
    we reject the sequence outright.

    A dictionary containing the following key-value pairs are returned:

    - "sequences": All proposed sequences.
    - "scores": All scores from the scoring function.
    - "accept": Whether the sequence was accepted as the new 'current sequence'
        on which new sequences are proposed.

    This can be turned into a pandas DataFrame.

    ### Parameters

    - `starter_sequence`: The starting sequence.
    - `n_steps`: Number of steps for the MC chain to walk.
    - `scoring_func`: Scoring function for a new sequence.
        It should only accept a string `sequence`.
    - `is_accepted_kwargs`: Dictionary of kwargs to pass into
        `is_accepted` function.
        See `is_accepted` docstring for more details.
    - `trust_radius`: Maximum allowed number of mutations away from
        starter sequence.
    - `propose_kwargs`: Dictionary of kwargs to pass into
        `propose` function.
        See `propose` docstring for more details.
    - `verbose`: Whether or not to print iteration number
        and associated sequence + score. Defaults to False

    ### Returns

    A dictionary with `sequences`, `accept` and `score` as keys.
    """
    current_sequence = starter_sequence
    current_score = scoring_func(sequence=starter_sequence)

    chain_data = defaultdict(list)
    chain_data["sequences"].append(current_sequence)
    chain_data["scores"].append(current_score)
    chain_data["accept"].append(True)

    for i in tqdm(range(n_steps)):
        new_sequence = propose(current_sequence, **propose_kwargs)
        new_score = scoring_func(sequence=new_sequence)

        default_is_accepted_kwargs = {"temperature": 0.1}
        default_is_accepted_kwargs.update(is_accepted_kwargs)
        accept = is_accepted(
            best=current_score,
            candidate=new_score,
            **default_is_accepted_kwargs,
        )

        # Check hamming distance
        if hamming_distance(starter_sequence, new_sequence) > trust_radius:
            accept = False

        # Determine acceptance
        if accept:
            current_sequence = new_sequence
            current_score = new_score

        # Record data.
        chain_data["sequences"].append(new_sequence)
        chain_data["scores"].append(new_score)
        chain_data["accept"].append(accept)
    chain_data["scores"] = np.hstack(chain_data["scores"])
    return chain_data

# Perform sampling
sampled_sequences = sample_one_chain(
    starting_sequence, n_steps=steps, scoring_func=scoring_func,
    trust_radius=trust, is_accepted_kwargs={'temperature': temperature}
)

# Extract the scores
scores = sampled_sequences.pop('scores')

# Create data frame
sampled_seqs_df = pd.DataFrame(sampled_sequences)

# Add score columns
sampled_seqs_df = pd.concat([
     sampled_seqs_df,
    pd.DataFrame({"Scores": scores})],
     axis=1
)

# Re-order columns
cols = sampled_seqs_df.columns.values
cols = list(cols[cols != 'accept']) + ['accept']

sampled_seqs_df = sampled_seqs_df[cols]

# Add step numbers
sampled_seqs_df['step'] = list(range(0, steps + 1))

# Save sampled sequences
sampled_seqs_df.to_csv(outfile, '\t', index=False)
