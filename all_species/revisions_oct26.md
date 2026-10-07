## Addressing reviewer comments on statistics

#### 1. AMOVA p-values
The reviewer noted that the "reported AMOVA / outlier p-values appear inverted or mis-specified" 
I found the issue in the code, which is that the AMOVA test reports 3 p-values, and my code was taking the first one. The first value was testing heterozygosity within individuals, but we needed to report the differences between populations. 
So, I replaced `pvalue[1]` in the code with `pvalue[length(pvalue)]`. 

```R
# amova FST
genI <- gl2gi(gl_bir)
pops_gi <- genI$pop; strata(genI) <- data.frame(pops = pops_gi)
p.amova <- poppr.amova(genI, ~pops)
amova.pvalues <- ade4::randtest(p.amova, nrepet = 9999)
overall_fst <- p.amova$statphi$Phi[length(p.amova$statphi$Phi)] # 0.003964456
# updated 10/7/26
print(amova.pvalues) # take row 3 from here for variation between pops 
overall_pval <- amova.pvalues$pvalue[length(amova.pvalues$pvalue)]
overall_pval # 9e-04

# fst but excluding BYC (n=1)
pop_sizes <- table(pop(gl_bir))
valid_pops <- names(pop_sizes[pop_sizes >= 2])
gl_valid <- gl_bir[pop(gl_bir) %in% valid_pops, ]
fst_result <- gl.fst.pop(gl_valid, nboots=1000, percent=95, nclusters=1)
fst_matrix <- fst_result$Fsts
fst_matrix[is.na(fst_matrix)] <- t(fst_matrix)[is.na(fst_matrix)]

# updated 10/7/26
pval_matrix <- fst_result$Pvalues
pval_matrix[is.na(pval_matrix)] <- t(pval_matrix)[is.na(pval_matrix)]
print("pairwise fst")
print(round(fst_matrix, 4))
print(pval_matrix)
## ok so noting the revision here, the results show significant comparisons between india-peru and india-mexico but not peru-mexico
# FST India vs Peru = 0.0067, p = 0 significant
# FST India vs Mexico = 0.0055, p = 0 significant 
# FST Peru vs Mexico = 0.0006, p= 0.278 not significant
```


#### 2. Outlier FST not tested
The reviewer said "the outlier FST values are not significance-tested"
They were significance tested, but the AMOVA had the same code issue as above, so I had to do the same fix.

```R
## outlier fst
table(ploidy(gl_outliers)) # 18 of the 22 inds are assigned ploidy 1 because there are no heterozygous calls across the 4 SNPs
ploidy(gl_outliers) <- 2 # set ploidy to 2 for all
nLoc(gl_outliers) # 4 outliers
pop(gl_outliers) <- popmap$pop[match(indNames(gl_outliers), popmap$sample)]
genI_out <- gl2gi(gl_outliers)
pops_gi_out <- genI_out$pop
strata(genI_out) <- data.frame(pops = pops_gi_out)

p.amova_out      <- poppr.amova(genI_out, ~pops)
amova.pvalues_out <- ade4::randtest(p.amova_out, nrepet = 9999)
overall_fst_out  <- p.amova_out$statphi$Phi[length(p.amova_out$statphi$Phi)]
# updated 10/7/26
print(amova.pvalues_out)
overall_pval_out <- amova.pvalues_out$pvalue[length(amova.pvalues_out$pvalue)]

# outlier amova fst and pval
round(overall_fst_out, 4) # fst = 0.3249
round(overall_pval_out, 4) # 1e-04
```

#### 3. Outlier bonferroni correction
The reviewer wanted more stringent significance testing for the outlier detection. 
Originally, we did add a q-value (FDR) correction at q < 0.01, but I also added a Bonferroni correction to report the number of loci. 

```R
# pcadapt for outlier loci
geno_bir <- read.pcadapt("/Users/madelineeppley/Desktop/manta26pub/birostris_22samp.vcf", type = "vcf")
obj <- pcadapt(geno_bir, K = 1)
pvals <- obj$pvalues
qvals <- qvalue(pvals)$qvalues
alpha <- 0.01
outliers <- which(qvals < alpha)
n_outliers <- length(outliers)
n_outliers # 4 outliers
outliers # 1081 1676 4995 8161

# updated 10/7/26
padj <- p.adjust(pvals, method = "bonferroni")
alphab <- 0.1
outliers_bonf <- which(padj < alphab)
length(outliers_bonf) #4 outliers
outliers # 1081 1676 4995 8161 same as above
```

#### 4. Pairwise p-values
The pairwise FST values originally came from `gl.fst.pop` but the p-values were also coming from the AMOVA test which had the same code issue as above. Instead of fixing the AMOVA, I'm just going to get rid of that here and pull the p-values from the `gl.fst.pop` output. 

```R
# pairwise FST excluding BYC (n=1)
pop_sizes_out  <- table(pop(gl_outliers))
valid_pops_out <- names(pop_sizes_out[pop_sizes_out >= 2])
gl_valid_out   <- gl_outliers[pop(gl_outliers) %in% valid_pops_out, ]

fst_result_out <- gl.fst.pop(gl_valid_out, nboots=1000, percent=95, nclusters=1)
fst_matrix_out <- fst_result_out$Fsts
fst_matrix_out[is.na(fst_matrix_out)] <- t(fst_matrix_out)[is.na(fst_matrix_out)]

# updated 10/7/26
pval_matrix_out <- fst_result_out$Pvalues
pval_matrix_out[is.na(pval_matrix_out)] <- t(pval_matrix_out)[is.na(pval_matrix_out)]
print("pairwise fst outlier snps")
print(round(fst_matrix_out, 4))
print(pval_matrix_out)
```

#### 5. Neutral loci FST 
There was an outdated set of neutral loci in this calculation, so I updated the script to use the correct set of neutral loci and re-calcluated with a fixed AMOVA. 

```R
# updated 10/7/26
stopifnot(nLoc(gl_bir) == nrow(vcf_bir@fix))
neutral <- setdiff(seq_len(nLoc(gl_bir)), outliers)
gl_neutral <- gl_bir[, neutral]
nLoc(gl_neutral) #12041
genI_neu <- gl2gi(gl_neutral)
strata(genI_neu) <- data.frame(pops = genI_neu$pop)
p.amova_neu <- poppr.amova(genI_neu, ~pops)
amova.pvalues_neu <- ade4::randtest(p.amova_neu, nrepet = 9999)
print(amova.pvalues_neu)
overall_fst_neu <- p.amova_neu$statphi$Phi[length(p.amova_neu$statphi$Phi)]
overall_pval_neu <- amova.pvalues_neu$pvalue[length(amova.pvalues_neu$pvalue)]
round(overall_fst_neu, 4)
round(overall_pval_neu, 4)
```

### Birostris
1. AMOVA p-values: Changed from p = 0.948 (not sig) to 0.0009 (sig).
2. Outlier FST: Changged from 0.368 to 0.00001 - this makes sense since the outliers SHOULD be sig different
3. Bonferroni correction: It was actually the same for birostris (n=4, same loci) that were outliers from both significance tests.
4. Pairwise p-values: India vs Peru and India vs Mexico are significant (p < 0.001) but Peru vs Mexico is not (p = 0.278)
5. Neutral loci FST: Updated to 12,041 neutral loci, FST changes to 0.0038, p = 0.0012 (sig)
