# Musical Babies // Ecological momentary assessment reveals causal effects of music enrichment on infant mood

This repo contains all files necessary for reproducing the article "Ecological momentary assessment reveals causal effects of music enrichment on infant mood" by Cho & Yurdum et al. 

**This repository is currently under construction. We will remove this warning when the repository is fully populated.**

## Description

You can find the following here:

- an R Markdown file that generates the manuscript
- data, analysis, and visualisation code to produce the results reported in the manuscript
- intervention materials 
- supplementary data and materials

**For assistance, please contact the corresponding authors: Eun Cho ([eun.cho@yale.edu](mailto:eun.cho@yale.edu)), Lidya Yurdum ([lidya.yurdum@yale.edu](mailto:lidya.yurdum@yale.edu)) and Samuel Mehr ([mehr@hey.com](mailto:mehr@hey.com)).**

## Anatomy of the repo

To render the paper, run the code in `/writing/manuscript.Rmd`. Running this script will generate the manuscript from saved outputs of analysis and visualization scripts. To run everything from the raw data (which can take a while to process), first start by running `/analysis/builder.R`, then `analysis.R`.

## Data and analysis code

All raw data files are in `/data`. Identifiable information have been removed.

Scripts for preprocessing the data are in `/analysis`.

Preprocessed data, interim datasets and the like are in `/results`.

### Visualisations

Visualisation code is in `/viz`, along with images and static data used for non-dynamic visualisations.

### Materials

Intervention materials are in `/materials`, and include video karaoke files provided to parents to aid in the singing intervention.
