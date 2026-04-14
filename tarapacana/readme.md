M. tarapacana
population structure results


<img width="8400" height="10800" alt="tarapacana_FINALFIGURE" src="https://github.com/user-attachments/assets/698d8374-2536-4bb9-aa83-daa226559068" />



## final analysis pipeline for bycatch (no geographic loc) samples
1. `01_tarapacana_param_opt.sh` parameter optimization script with m=3 fixed, tested M=2,3,4,5, and n-M-1/M/M+1. Removed BYC_RMB_56 first because it is closely related to 57. Used a single population map called opt to avoid -p filtering issues. Pipline was ustacks -> cstacks -> sstacks -> tsv2bam -> gstacks -> populations. output went into /opt directories.
2. `02_tarapacana_continue.sh` re-ran interrupted optimization downstream steps (previous job timed out). Looped through existing /opt directories and picked up with running sstacks and downstream steps if needed.
3. `03_tarapacana_final_run.sh` final pipeline with the optimal params determed by the opt script. This was (m=3, M=2, n=2). Included 14 individuals. The populations module was originally run with the tarapacana_pop_map_3, which had 4 populations (EST, OFF, CTR, WST), but this is replaced by re-running populations in the following script. This used the same full pipeline as above and used the VCFtools filtering with minDP 10, max-missing 0.8, and removed individuals with >40% data missing. the output was in the final_m3_M2_n2/ directory.
4. `04_run_populations_gal.sh` re-runs just the populations module and the VCF tools filtering with the new groupings. input was the existing STACKS catalog from the final_m3_M2_n2/ folder. The popmap_gal.txt with new assignments (EAST=5, GAL=4, CTR=3, WST=1). n=13 samples remained after filtering (BYC_RMT_59 was removed during VCFtools filtering). This is the final dataset.
