<img src="https://github.com/AnantharamanLab/METABOLIC/blob/master/METABOLIC.jpg" width="85%">

# Note:

The `dev` branch is where I'm currently converting this pipeline into a Nextflow workflow.


# METABOLIC
**MET**abolic **A**nd **B**ioge**O**chemistry ana**L**yses **I**n mi**C**robes  

This software enables the prediction of metabolic and biogeochemical functional trait profiles to any given genome datasets. These genome datasets can either be metagenome-assembled genomes (MAGs), single-cell amplified genomes (SAGs) or isolated strain sequenced genomes. 

# Purpose
This pipeline primary purpose is to summarize the putative involvement of genomes-level omics data in biogeochemical pathways.
Two primary options are available: the `genome` screens for the presence or absense of genes in specific pathways, and summarizes the information into tables and figures for the C, N and S cycles.

The `community` option is the same as the `genome` option, but furthers add the requirement to provide a set of sequence reads in order to quantity the abundance of the genomes or the transcript-level expression (in the case of providing transcriptional reads). The `community` version of the pipeline creates additional figures, such as the functional network diagram, an alluvial plot showing which taxa are involved in which cycle, and a figure showing potential metabolic hand-offs. In this version, all the figures are generated at the `community` level for an overview, rather than at the `genome` level. 

# Updates to the workflow
In this updated version of the pipeline, the workflow has been converted into a `Nextflow` pipeline to offer greater reproducibility and scalability.

To run this pipeline, you will need to:
- Obtain a copy of this github repository
- Set up the proper databases:
    - HMM profiles
    - KoFAM
    - MEROPS [Nextflow manages this]
    - DBCAN2 (Cazyme) [Nextflow manages this]
    - EGGNOG-Mapper (*new!*) [Nextflow manages this]
- Create a project where to run the Nextflow program

# Dependencies
You will need to have `nextflow` and `conda` installed on your system. Each step of the workflow is broken down into individual processes under `modules/local/process.nf`, and managed via `conda`. 

Nextflow documentation: https://www.nextflow.io/docs/latest/install.html
Miniconda3 documentation: https://www.anaconda.com/docs/getting-started/miniconda/main

# Using Nextflow
Any parameter in `nextflow.config` can be given as a command line option using double dashes (`--parameter`). For example, if you already have your databases set-up and want to reuse them, using `--METABOLIC_db_path /path/to/METABOLIC_db`. 
We also recommend using the `-resume` option (no double dash), to resume runs. For example, if you have ran most steps in the workflow, stopped the run and want to re-run this again, nextflow will already know which steps are completed (cached), and restart at the appropriate point.

# Logs
A runtime log and general report will be created in your project folder, and named `report.html` and `timeline.html`. This Nextflow-style report can be opened with any web browser, and give you an idea of which steps took longest to run, how long they ran for, which steps succeeded vs. failed, etc.

# Detailed documentation

[TBD]

# Citation

Previous versions:

If you are using this program, please consider citing our paper, available at [Microbiome](https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-021-01213-8):
```
Zhou, Z., Tran, P.Q., Breister, A.M. et al. METABOLIC: high-throughput profiling of microbial genomes for functional traits, metabolism, biogeochemistry, and community-scale functional networks. Microbiome 10, 33 (2022). https://doi.org/10.1186/s40168-021-01213-8
```

