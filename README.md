## Summary

This repository contains all code and figures for the manuscript: 

"Genetic population structure and conservation implications for manta and devil rays in the Eastern Pacific Ocean"

Authors: Melissa R. Cronin, Madeline G. Eppley, Remy Gatins, Giacomo Bernardi, Maria de Lourdes Torres, Daniel Fernando, Michel Guererro, Jon Lopez, Gala Moreno, Marta D. Palacios, Alessandro Ponzo, Salvador Siu, Joshua D. Stewart, Daniel B. Wright, Kelly M. Zilliacus, Donald A. Croll

<img width="4000" height="2250" alt="Presentation3" src="https://github.com/user-attachments/assets/caa9044c-7fff-4d03-9584-1fb7567b7264" />

## Structure
#### STACKS pipeline
`01_species_param_opt.sh` parameter optimization script with m=3 fixed, tested M=2,3,4,5, and n-M-1/M/M+1. Used a single population map called opt to avoid -p filtering issues. Pipline was ustacks -> cstacks -> sstacks -> tsv2bam -> gstacks -> populations. Output into /opt directories.

`02_species_continue.sh` *optional* for if parameterization pipeline was interrupted during optimization downstream steps (e.g., previous job timed out). Looped through existing /opt directories and picked up with running sstacks and downstream steps if needed.

`03_species_final_run.sh` final pipeline with the optimal STACKS params determed by the opt script (e.g., m=3, M=2, n=2). Includes vcftools filtering with minDP 10, max-missing 0.8, and removed individuals with >40% data missing. the output was in the optimal param directory per species (e.g., `final_m3_M2_n2/`) directory.



## Data
Upon final acceptance, all sequence data will be archived on the NCBI sequence archive and the project accession numbers will be listed here and in the published work. 
