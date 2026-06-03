## pcadapt with qvalue/FDR correction

I realized that I missed a line of code and currently our outliers are just using raw p-vals from pcadapt. Here's how I corrected the code:

```R
# pcadapt 
geno_tar <- read.pcadapt("/Users/madelineeppley/Desktop/manta/finalvcfs/tarapacana_final.vcf", type="vcf")
obj <- pcadapt(geno_tar, K=1)
pvals <- obj$pvalues
qvals <- qvalue(pvals)$qvalues
alpha <- 0.01
outliers <- which(qvals < alpha)
n_outliers <- length(outliers)
n_outliers
```
To get structure plots with pop labels quickly
```
> structure_outlier
> layer3 <- wrap_plots(structure_outlier, pop_label_plot, ncol=1, heights=c(1,0.3))
> layer3
```

Re-running pcadapt yields the following numbers of q-value/FDR corrected outliers:

Number of outliers with FDR/qval correction
- Birostris: 4 (originally 165)
- Munkiana: 120 (originally 211)
- Thurstoni: 13 (originally 41)
- Mobular: 260 (originally 1,003)
- Tarapacana: 116 (originally 1,922)


## New mobular plot
This looks extremely similar to the current outlier PCA plot
<img width="4800" height="3600" alt="mobular_PCA_outlierSNPs" src="https://github.com/user-attachments/assets/94fe714c-2407-4f24-b10c-87bc3407038e" />

<img width="1370" height="602" alt="Screenshot 2026-06-03 at 11 51 26 AM" src="https://github.com/user-attachments/assets/6af92a7b-bbb5-4f11-9dcf-f869f2f34c3d" />


## New munkiana plot
Same vibes essentially 
<img width="4800" height="3600" alt="munkiana_PCA_outlierSNPs" src="https://github.com/user-attachments/assets/dd7cb320-484f-470e-8746-9246fc434f54" />

<img width="889" height="374" alt="Screenshot 2026-06-03 at 11 53 30 AM" src="https://github.com/user-attachments/assets/9ff2f17c-da19-4a45-9ab6-5039a55cdc41" />


## New tarapacana plot
Maybe we see a little more structure here with the Western sample? I think this could aid our discussion
<img width="4800" height="3600" alt="tarapacana_PCA_outlierSNPs" src="https://github.com/user-attachments/assets/5bb8fb78-8bb0-4082-9e8c-b43030faf2a7" />


## New birostris plots
Note only 4 outliers here.

<img width="4800" height="3600" alt="birostris_PCA_outlierSNPs" src="https://github.com/user-attachments/assets/d4eab934-79a6-4996-a1ae-b0710314dfb0" />

<img width="1586" height="657" alt="Screenshot 2026-06-03 at 11 45 28 AM" src="https://github.com/user-attachments/assets/d8677b6e-ac2a-4d18-900f-53fe624e0986" />

