# Protein-Design

### Description

This repo gives you the workflow for designing new alpha-amylase variants with improved activity.
We use the structure of alpha-amylase from PDB iD [1BAG](https://www.rcsb.org/structure/1BAG) and 
use [ESM3](https://github.com/evolutionaryscale/esm) to generate candidate designs for further evaluation.
We also build an ML model to predict alpha-amylase activity from the amino acid sequences to score our designs.

The overall workflow is as below
<p align="center">
<img src="Data/figs/workflow.jpg" />
</p>

----

### Details about the Repo

#### Initial Analysis

A preliminary analysis was conducted on the given dataset (Data/sequences.csv). 
Some sequences were found to be different from alpha-amylases so they were removed.

[notebook: load_and_analyse.ipynb](src/notebooks/load_analyse.ipynb) - we create fasta files and eliminate suspicious sequences.

[notebook: validate_input_sequences.ipynb](src/notebooks/validate_input_sequences.ipynb) - here we see if the suspicious sequences are actually not alpha-amylases. We do this using hmm profiles for the domains in alpha-amylases and using hmmsearch on the suspected sequences.

