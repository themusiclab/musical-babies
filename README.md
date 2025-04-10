# Ecological momentary assessment reveals causal effects of music enrichment on infant mood

This repo contains all files necessary for reproducing the article "Ecological momentary assessment reveals causal effects of music enrichment on infant mood".

## Description

You can find the following here:

- an R Markdown file that generates the manuscript
- data, analysis, and visualisation code to produce the results reported in the manuscript
- intervention materials 
- supplementary data and materials


To set the working directory, open `musical-babies.Rproj`. All analysis scripts should be run from within this project.

## Anatomy of the repo

To render the paper, run the code in `/writing/MIPH_childdev.Rmd`. Running this script will generate the manuscript from saved outputs of the analysis script. To run everything from the cleaned data, run `analysis.R`.

## Data and analysis code

All data files are in `/data`. Identifiable information have been removed.

Analysis scripts are in `/analysis`.

Analysis results are in `/results`.

### Visualisations

Images used for non-dynamic visualisations are in `/viz`.

[![DOI](https://zenodo.org/badge/787304983.svg)](https://doi.org/10.5281/zenodo.15181606)

### Materials

Intervention materials (i.e., video karaoke files provided to parents to aid in the singing intervention) are in `/materials`.
