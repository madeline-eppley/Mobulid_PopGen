## workflow plan for each species
- STACKS, test parameter space and plot outcomes of # of mapped loci for all species, this filters for LD
- filter using vcftools for missingness at the individual and SNP level, MAF
- test for relatedness and remove highly related individuals
- filtered vcfs -> genlight objects, etc.
- report coverage stats
- input filtered vcfs into R -> PCA, ancestry K 1-5 show all plots, cross-entropy plots
- Fst heatmaps between sampling sites

other outputs
what kind of report do you want for outliers? 


For species with bycatch samples: 
- sometimes I think it will make sense to pre-group the sample with the closest coastal sampling group, but this could also lead to an effect where 
