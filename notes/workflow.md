## workflow plan
- STACKS, test parameter space and plot outcomes of # of mapped loci for all species
- filter using vcftools for missingness at the individual and SNP level
- test for relatedness and remove highly related individuals
- filter for LD in PLINK using -indep with a sliding window fo 50 SNPs (also could test 100?)
- filtered vcfs -> genlight objects, produce STRUCTURE inputs
- PCA, ancestry K 1-5 show all plots, cross-entropy plots
- 

## filtering with vcftools 
From the [Humble 2023 supplementary file](https://onlinelibrary.wiley.com/action/downloadSupplement?doi=10.1111%2Fmec.17220&file=mec17220-sup-0001-AppendixS1.pdf): VCFtools v1.13 in order to remove
low quality sites and individuals. For this, genotypes with a depth of coverage less than 6,
SNPs called in less than 40% individuals and individuals with more than 45% missing data
were removed. We then filtered SNPs more stringently by removing those called in less than
60% of individuals and with a depth of coverage less than 25 or greater than 376 (equivalent
to 95% quantiles for mean depth of coverage).

## filtering for LD
From the [Humble 2023 supplementary file](https://onlinelibrary.wiley.com/action/downloadSupplement?doi=10.1111%2Fmec.17220&file=mec17220-sup-0001-AppendixS1.pdf): "SNPs were pruned for linkage
disequilibrium using the -indep function in PLINK v1.9 with a sliding window of 50 SNPs, a
step size of 5 SNPs and variance inflation threshold of 2. These datasets were used for analysis of population structure and contemporary
migration."

## optimizing loci output from stacks 
figure shows the change in the number of polymorhpic loci present in at least 80% of the samples for increasing values of M and n assembly parameters in stacks. 
In Humble 2023 the number of assembled loci reached peak at M=3 and n=4. 

## relatedness
also ran the KING-robust kinship for empirical relatedness for all individaul pairwise comaprisons

## ancestry proportions for K=1 to K=8 were shown in a composite plot
- for birostris, no clear designations at K=2, cross-validation error increases after 1 with no knee
- for alfredi, clear designations at K=2 through K=5, cross-validation error shows knee in the plot at K=4

- ## Fst estimates for all pariwise population comparisons

