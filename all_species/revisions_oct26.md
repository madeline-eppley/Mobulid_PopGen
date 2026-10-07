## Addressing reviewer comments on statistics

#### AMOVA p-values
The reviewer noted that the "reported AMOVA / outlier p-values appear inverted or mis-specified" 
I found the issue in the code, which is that the AMOVA test reports 3 p-values, and my code was taking the first one. The first value was testing heterozygosity within individuals, but we needed to report the differences between populations. 
So, I replaced `pvalue[1]` in the code with `pvalue[length(pvalue)]`. For instance, this changed birostris from p = 0.948 (not sig) to 0.0009 (sig).


#### Outlier FST not tested
The reviewer said "the outlier FST values are not significance-tested"
They were significance tested, but the AMOVA had the same code issue as above, so I had to do the same fix. This makes the outlier FST value significant for birostris

#### Outlier bonferroni correction
The reviewer wanted more stringent significance testing for the outlier detection. 
Originally, we did add a q-value (FDR) correction at q < 0.01, but I also added a Bonferroni correction to report the number of loci. It was actually the same for birostris (n=4, same loci) that were outliers from both significance tests.

#### Pairwise p-values
The pairwise FST values originally came from `gl.fst.pop` but the p-values were also coming from the AMOVA test which had the same code issue as above. 
