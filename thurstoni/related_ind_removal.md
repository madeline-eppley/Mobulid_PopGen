## MGE notes 6/10/2026

We have two individuals that are closely related in the thurstoni dataset. The question is do we remove one for keep both for the publication figure. This depends on if the two individuals biasing population structure?

## ISSUES 
1. When we run the analysis with 1 highly related individual removed, we get 0 q-value 0.01 corrected outliers.
2. When we re-ran with the q-value of 0.01 cutoff and all individuals, we get 13 outliers. However, we can't calculate FST on these because the missing data exceeds 10% across all individuals for these SNPs. The standard cutoff value for missingness is 5% to calculate FST, so we're really out of range.
3. The thurstoni dataset is by far the smallest (1.1k total SNPs, 13 outliers with q-value 0.01 cutoff, higher levels of missingness across individuals). I relaxed a vcftools missingness parameter to increase our SNP yield for this species already. 

### K = 3 all SNPs when the 2 highly related individuals are kept
- I think this is probably the version that we should keep as the main figure of the manuscript
<img width="8400" height="10800" alt="thurstoni_FINALFIGURE" src="https://github.com/user-attachments/assets/a6973cb0-32cd-47a3-8f0a-854dc6823482" />


### Population structure results with 1 highly related individual removed:
- There's really not a huge change here in how we would interpret population structure. I could add just the PCA and the DAPC to the supplement, if we would like. 
<img width="8400" height="10800" alt="thurstoni_FINALFIGURE_v2" src="https://github.com/user-attachments/assets/c8e6647d-11df-49ff-95fe-ac1b8b4bb6d1" />




