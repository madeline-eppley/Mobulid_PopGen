## Addressing reviewer comments on statistics

#### 1. AMOVA p-values
The reviewer noted that the "reported AMOVA / outlier p-values appear inverted or mis-specified" 
I found the issue in the code, which is that the AMOVA test reports 3 p-values, and my code was taking the first one. The first value was testing heterozygosity within individuals, but we needed to report the differences between populations. 
So, I replaced `pvalue[1]` in the code with `pvalue[length(pvalue)]`. 


#### 2. Outlier FST not tested
The reviewer said "the outlier FST values are not significance-tested"
They were significance tested, but the AMOVA had the same code issue as above, so I had to do the same fix.

#### 3. Outlier bonferroni correction
The reviewer wanted more stringent significance testing for the outlier detection. 
Originally, we did add a q-value (FDR) correction at q < 0.01, but I also added a Bonferroni correction to report the number of loci. 

#### 4. Pairwise p-values
The pairwise FST values originally came from `gl.fst.pop` but the p-values were also coming from the AMOVA test which had the same code issue as above. Instead of fixing the AMOVA, I'm just going to get rid of that here and pull the p-values from the `gl.fst.pop` output. 

#### 5. Neutral loci FST 
There was an outdated set of neutral loci in this calculation, so I updated the script to use the correct set of neutral loci and re-calcluated with a fixed AMOVA. 

### Birostris
1. AMOVA p-values: Changed from p = 0.948 (not sig) to 0.0009 (sig).
2. Outlier FST: Changged from 0.368 to 0.00001 - this makes sense since the outliers SHOULD be sig different
3. Bonferroni correction: It was actually the same for birostris (n=4, same loci) that were outliers from both significance tests.
4. Pairwise p-values: India vs Peru and India vs Mexico are significant (p < 0.001) but Peru vs Mexico is not (p = 0.278)
5. Neutral loci FST: Updated to 12,041 neutral loci, FST changes to 0.0038, p = 0.0012 (sig)
