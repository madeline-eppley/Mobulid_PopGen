# Birostris RADseq analysis

In this dataset, we have 36 RAD-sequenced samples from India, Peru, Mexico, and one bycatch sample from the high seas.

Mobula birostris, or the giant oceanic manta ray, is the largest ray in the world :o 

It is an endangered species in tropical and temperate waters across the world. let's get into some genetics!!!

<img width="554" height="554" alt="image" src="https://github.com/user-attachments/assets/ec6d715f-909b-4505-8c5b-4aa676f1c872" />


### fix the pop map
move our birostris pop map file into the working dir with 
```bash
cp /projects/gatins/2025_Mobulid_UCSC/RAD_all_combined_bycatch/pop_map_birostris /projects/gatins/2025_Mobulid/birostris/pop_map_birostris
```

now clean up the missingness
```bash
(base) [eppley.m@explorer-02 RAD_all_combined_bycatch]$ cat pop_map_birostris 
BYC_RMB_01	BYC
IN_1_MB	IND
IN_18_MB_B	IND
IN_2_MB	IND
IN_4_MB	IND
IN_5_MB	IND
IN_6_MB	IND
IN_7_MB	IND
PER_001_MB	PER
PER_003_MB	PER
PER_004_MB	PER
PER_005_MB	PER
PER_006_MB	PER
PER_007_MB	PER
PER_008_MB	PER
PER_DZW81_4_MB	PER
REV_1_MB	REV
REV_10_MB	REV
REV_13_MB	REV
REV_14_MB_B	REV
REV_15_MB	REV
REV_17_MB	REV
REV_18_MB	REV
REV_19_MB	REV
REV_20_MB_B	REV
REV_5_MB	REV
```

cleaned version: 
```bash
BYC_RMB_01	BYC
IN_1_MB	IND	IND
IN_18_MB_B	IND
IN_2_MB	IND	IND
IN_4_MB	IND	IND
IN_5_MB	IND	IND	
IN_6_MB	IND	IND
IN_7_MB	IND	IND
PER_001_MB	PER
PER_003_MB	PER
PER_004_MB	PER
PER_005_MB	PER
PER_006_MB	PER
PER_007_MB	PER
PER_008_MB	PER
PER_DZW81_4_MB	PER
REV_1_MB	REV
REV_10_MB	REV
REV_13_MB	REV
REV_14_MB_B	REV
REV_15_MB	REV
REV_17_MB	REV
REV_18_MB	REV
REV_19_MB	REV
REV_20_MB_B	REV
REV_5_MB	REV
```

## the bycatch sample
when we run stacks, if I set the bycatch sample to be it's own population, it will compare the loci in that individual to all others and could have an effect on overall loci retention because it will only keep loci present in that individual. therefore for the intital run i'm going to look a the collection lat/long for this sample and assign it to the closest population if it makes reasonable sense. 

collection information for individaul BYC_RMB_01
Número	Fecha	Especie	Sexo	Latitud	Longitud	ID
1	10/17/19	RMB	H	0428S	09901W	BYC_RMB_01

decimal lat/long is : -4.466700, -99.016700

which falls a bit off the coast of Peru. Let's assign that individual to the PER population and go from there. 
<img width="1155" height="643" alt="Screenshot 2025-11-13 at 10 55 13 AM" src="https://github.com/user-attachments/assets/e33e2ac5-0e3d-4eb6-ad94-a731d3d4ce8c" />

updates pop map file ... which brings us to 3 popualtions, PER, IND, and REV.

### pop map file formatting
because i just edited the pop map by hand and pressed tab (lazy) i got an error about the formatting. I had created a new column with tab, but a pop map is limited to 2 columns. in case I need this again, here's how i fixed this issue: 

```bash
(base) [eppley.m@explorer-02 birostris]$ cat -A pop_map_birostris | head
BYC_RMB_01^IPER$
IN_1_MB^IIND^IIND$
IN_18_MB_B^IIND$
IN_2_MB^IIND^IIND$
IN_4_MB^IIND^IIND$
IN_5_MB^IIND^IIND^I$
IN_6_MB^IIND^IIND$
IN_7_MB^IIND^IIND$
PER_001_MB^IPER$
PER_003_MB^IPER$
(base) [eppley.m@explorer-02 birostris]$ awk '{print $1 "\t" $2}' pop_map_birostris > pop_map_birostris.clean && mv pop_map_birostris.clean pop_map_birostris
(base) [eppley.m@explorer-02 birostris]$ cat -A pop_map_birostris | head
BYC_RMB_01^IPER$
IN_1_MB^IIND$
IN_18_MB_B^IIND$
IN_2_MB^IIND$
IN_4_MB^IIND$
IN_5_MB^IIND$
IN_6_MB^IIND$
IN_7_MB^IIND$
PER_001_MB^IPER$
PER_003_MB^IPER$
(base) [eppley.m@explorer-02 birostris]$ awk '{print NF}' pop_map_birostris | sort | uniq -c
     26 2
```

at the end there, we have 26 samples x 2 cols! sweet


## full pipeline
now let's execute our full stacks pipeline! we are using our standard settings here that we established after running the tarapacana dataset. 

```bash
#!/bin/bash

#SBATCH --partition=lotterhos
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --time=48:00:00
#SBATCH --job-name=stacks_birostris
#SBATCH --output=stacks_birostris_%j.log
#SBATCH --error=stacks_birostris_%j.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=eppley.m@northeastern.edu


export PATH=/projects/gatins/programs_explorer/stacks_2.68/bin:$PATH


echo "Sample list:"
cat /projects/gatins/2025_Mobulid/birostris/pop_map_birostris

echo "Starting denovo_map.pl"

denovo_map.pl \
  -m 3 \
  -M 3 \
  -n 2 \
  -T 32 \
  -o /projects/gatins/2025_Mobulid/birostris \
  --popmap /projects/gatins/2025_Mobulid/birostris/pop_map_birostris \
  --samples /projects/gatins/2025_Mobulid_UCSC/RAD_all_combined_bycatch/trimmed90 \
  -X "populations:-r 0.8 -p 3 --min-maf 0.05 --write-single-snp --vcf --genepop --structure --fstats --hwe -t 30"


ls -l /projects/gatins/2025_Mobulid/birostris

echo "populations module done"

echo "starting vcftools"
module load vcftools # this is just globally available on explorer, slay

# input VCF from stacks population run
INPUT_VCF="/projects/gatins/2025_Mobulid/birostris/populations.snps.vcf"

OUTDIR="/projects/gatins/2025_Mobulid/birostris"
mkdir -p ${OUTDIR}

# filter by minimum depth per genotype (minDP = 10)
vcftools --vcf ${INPUT_VCF} \
         --minDP 10 \
         --recode --recode-INFO-all \
         --out ${OUTDIR}/minDP10

# filter for sites present in >= 80% of individuals
vcftools --vcf ${OUTDIR}/minDP10.recode.vcf \
         --max-missing 0.8 \
         --recode --recode-INFO-all \
         --out ${OUTDIR}/minDP10_maxmiss0.8

# remove individuals with >40% missing data
# 1: compute missingness per individual
vcftools --vcf ${OUTDIR}/minDP10_maxmiss0.8.recode.vcf \
         --missing-indv \
         --out ${OUTDIR}/missingness

# 2: generate a list of individuals to remove
awk '$5 > 0.4 {print $1}' ${OUTDIR}/missingness.imiss > ${OUTDIR}/remove_individuals.txt

# now filter out individuals
vcftools --vcf ${OUTDIR}/minDP10_maxmiss0.8.recode.vcf \
         --remove ${OUTDIR}/remove_individuals.txt \
         --recode --recode-INFO-all \
         --out ${OUTDIR}/minDP10_maxmiss0.8_filtInd
```

## results
tldr; we are keeping 22 individuals and 11.8k SNPs as the final dataset, here are the results: 

#### gstacks output
```bash
Genotyped 178705 loci:
  effective per-sample coverage: mean=20.3x, stdev=11.1x, min=7.4x, max=47.6x
  mean number of sites per locus: 90.0
  a consistent phasing was found for 50638 of out 62583 (80.9%) diploid loci needing phasing
```

```
Depths of Coverage for Processed Samples:
BYC_RMB_01; loci assembled: 130803; depth: 36.83x; max: 811x; number of reads: 4782079 (85.4%)
IN_1_MB; loci assembled: 121520; depth: 21.11x; max: 490x; number of reads: 2549088 (82.3%)
IN_18_MB_B; loci assembled: 81343; depth: 8.64x; max: 314x; number of reads: 698981 (79.3%)
IN_2_MB; loci assembled: 101626; depth: 12.05x; max: 314x; number of reads: 1216263 (81.0%)
IN_4_MB; loci assembled: 108573; depth: 10.47x; max: 398x; number of reads: 1130003 (81.7%)
IN_5_MB; loci assembled: 123697; depth: 15.48x; max: 367x; number of reads: 1903358 (83.6%)
IN_6_MB; loci assembled: 126167; depth: 17.20x; max: 395x; number of reads: 2156576 (81.4%)
IN_7_MB; loci assembled: 124105; depth: 18.53x; max: 479x; number of reads: 2285649 (85.0%)
PER_001_MB; loci assembled: 131413; depth: 39.71x; max: 942x; number of reads: 5182634 (84.1%)
PER_003_MB; loci assembled: 128513; depth: 26.96x; max: 589x; number of reads: 3443386 (85.1%)
PER_004_MB; loci assembled: 124339; depth: 18.73x; max: 483x; number of reads: 2314053 (84.4%)
PER_005_MB; loci assembled: 129786; depth: 38.94x; max: 839x; number of reads: 5018910 (84.5%)
PER_006_MB; loci assembled: 129702; depth: 23.49x; max: 507x; number of reads: 3028279 (83.7%)
PER_007_MB; loci assembled: 129146; depth: 21.03x; max: 643x; number of reads: 2698557 (84.3%)
PER_008_MB; loci assembled: 108916; depth: 9.24x; max: 306x; number of reads: 1001335 (83.1%)
PER_DZW81_4_MB; loci assembled: 129910; depth: 45.07x; max: 1335x; number of reads: 5814942 (84.3%)
REV_1_MB; loci assembled: 97320; depth: 7.60x; max: 304x; number of reads: 735713 (80.7%)
REV_10_MB; loci assembled: 112843; depth: 9.52x; max: 204x; number of reads: 1067764 (80.9%)
REV_13_MB; loci assembled: 130746; depth: 22.90x; max: 407x; number of reads: 2975946 (82.2%)
REV_14_MB_B; loci assembled: 123354; depth: 16.90x; max: 490x; number of reads: 2072218 (85.4%)
REV_15_MB; loci assembled: 119381; depth: 13.36x; max: 408x; number of reads: 1585774 (83.8%)
REV_17_MB; loci assembled: 121974; depth: 15.94x; max: 475x; number of reads: 1933145 (84.5%)
REV_18_MB; loci assembled: 118268; depth: 10.67x; max: 307x; number of reads: 1254545 (83.3%)
REV_19_MB; loci assembled: 129404; depth: 28.21x; max: 630x; number of reads: 3626545 (84.2%)
REV_20_MB_B; loci assembled: 122286; depth: 13.96x; max: 374x; number of reads: 1697638 (85.8%)
REV_5_MB; loci assembled: 96905; depth: 7.30x; max: 266x; number of reads: 704001 (79.0%)
```

#### after filtering for minDP 10
```
After filtering, kept 26 out of 26 Individuals
Outputting VCF file...
After filtering, kept 22958 out of a possible 22958 Sites
Run Time = 1.00 seconds
```

#### after filtering for max-missing 0.8
```
After filtering, kept 26 out of 26 Individuals
Outputting VCF file...
After filtering, kept 11843 out of a possible 22958 Sites
Run Time = 1.00 seconds
```

#### after filtering for missingness at the individual level 
```bash
After filtering, kept 22 out of 26 Individuals
Outputting VCF file...
After filtering, kept 11843 out of a possible 11843 Sites
```

### removed individuals
luckily, we have an even spread of individuals across populations that didn't pass our quality filtering. these will be removed from the analysis.
```
(base) [eppley.m@explorer-02 birostris]$ cat remove_individuals.txt 
INDV
IN_18_MB_B
PER_008_MB
REV_1_MB
REV_5_MB
```

### transferring vcf to local R
scp eppley.m@login.explorer.northeastern.edu:/projects/gatins/2025_Mobulid/birostris/minDP10_maxmiss0.8_filtInd.recode.vcf ~/Desktop/manta/birostris


### Initial data exploration in R 
### PCA 
I think I could see some regional groupings here ... need to import pops into R and plot PCA with those
Notably our bycatch sample caught off of peru looks like it's grouping near other samples from Peru! super cool :D

<img width="1413" height="1119" alt="image" src="https://github.com/user-attachments/assets/d3898412-3dd2-4310-ae90-a941ed484f85" />

<img width="1413" height="1119" alt="image" src="https://github.com/user-attachments/assets/dd044ea8-0ed5-4747-972b-9ba87dbe3b74" />

<img width="1413" height="1119" alt="image" src="https://github.com/user-attachments/assets/c1d2e8e8-3ecb-49aa-9392-a8114307fc3c" />

### relatedness
it kind of looks like REV_10 is more closely related to everything else ....??

<img width="1980" height="1494" alt="image" src="https://github.com/user-attachments/assets/e46aba31-30f9-42ab-bf6c-0adc2b255105" />


### manhattan
179 outliers
<img width="1980" height="1494" alt="image" src="https://github.com/user-attachments/assets/a86e553e-e4d1-4741-b035-ada8ae5c9307" />



## improving visualizations in R
i think that our workflow looks reasonably good at this point and we can proceed with our final visulization goals for birostris!

step 1:  get the pop map file into R
`scp eppley.m@login.explorer.northeastern.edu:/projects/gatins/2025_Mobulid/birostris/pop_map_birostris ~/Desktop/manta/birostris`

here's our workflow:
1. filtered vcfs -> genlight objects
2. PCA
3. LEA ancestry analysis, and ancestry K 1-5 show all plots plus the cross-entropy plot
4. Fst heatmaps between sampling sites



# running the byc sample as it's own population
Remy and I met and she had the idea to run the bycatch sample as it's own population to compare the results when we group it with the Peru pop (as done above). let's make a copy of our current pop map, switch analysis over to a new folder, and then change the assigned population of our bycatch sample from PER to BYC.

```bash
cp /projects/gatins/2025_Mobulid/birostris/pop_map_birostris /projects/gatins/2025_Mobulid/birostris/byc_own_pop/pop_map_birostris
```

```
(base) [eppley.m@explorer-01 byc_own_pop]$ cat pop_map_birostris 
BYC_RMB_01	BYC
IN_1_MB	IND
IN_18_MB_B	IND
IN_2_MB	IND
IN_4_MB	IND
IN_5_MB	IND
IN_6_MB	IND
IN_7_MB	IND
PER_001_MB	PER
PER_003_MB	PER
PER_004_MB	PER
PER_005_MB	PER
PER_006_MB	PER
PER_007_MB	PER
PER_008_MB	PER
PER_DZW81_4_MB	PER
REV_1_MB	REV
REV_10_MB	REV
REV_13_MB	REV
REV_14_MB_B	REV
REV_15_MB	REV
REV_17_MB	REV
REV_18_MB	REV
REV_19_MB	REV
REV_20_MB_B	REV
REV_5_MB	REV
```


### ok now let's rework the pipeline to work out of our new directory and use the new pop map

```bash
#!/bin/bash

#SBATCH --partition=lotterhos
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --time=48:00:00
#SBATCH --job-name=stacks_birostris
#SBATCH --output=stacks_birostris_%j.log
#SBATCH --error=stacks_birostris_%j.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=eppley.m@northeastern.edu


export PATH=/projects/gatins/programs_explorer/stacks_2.68/bin:$PATH


echo "Sample list:"
cat /projects/gatins/2025_Mobulid/birostris/byc_own_pop/pop_map_birostris

echo "Starting denovo_map.pl"

denovo_map.pl \
  -m 3 \
  -M 3 \
  -n 2 \
  -T 32 \
  -o /projects/gatins/2025_Mobulid/birostris/byc_own_pop \
  --popmap /projects/gatins/2025_Mobulid/birostris/byc_own_pop/pop_map_birostris \
  --samples /projects/gatins/2025_Mobulid_UCSC/RAD_all_combined_bycatch/trimmed90 \
  -X "populations:-r 0.8 -p 3 --min-maf 0.05 --write-single-snp --vcf --genepop --structure --fstats --hwe -t 30"


ls -l /projects/gatins/2025_Mobulid/birostris/byc_own_pop

echo "populations module done"

echo "starting vcftools"
module load vcftools # this is just globally available on explorer, slay

# input VCF from stacks population run
INPUT_VCF="/projects/gatins/2025_Mobulid/birostris/byc_own_pop/populations.snps.vcf"

OUTDIR="/projects/gatins/2025_Mobulid/birostris/byc_own_pop"
mkdir -p ${OUTDIR}

# filter by minimum depth per genotype (minDP = 10)
vcftools --vcf ${INPUT_VCF} \
         --minDP 10 \
         --recode --recode-INFO-all \
         --out ${OUTDIR}/minDP10

# filter for sites present in >= 80% of individuals
vcftools --vcf ${OUTDIR}/minDP10.recode.vcf \
         --max-missing 0.8 \
         --recode --recode-INFO-all \
         --out ${OUTDIR}/minDP10_maxmiss0.8

# remove individuals with >40% missing data
# 1: compute missingness per individual
vcftools --vcf ${OUTDIR}/minDP10_maxmiss0.8.recode.vcf \
         --missing-indv \
         --out ${OUTDIR}/missingness

# 2: generate a list of individuals to remove
awk '$5 > 0.4 {print $1}' ${OUTDIR}/missingness.imiss > ${OUTDIR}/remove_individuals.txt

# now filter out individuals
vcftools --vcf ${OUTDIR}/minDP10_maxmiss0.8.recode.vcf \
         --remove ${OUTDIR}/remove_individuals.txt \
         --recode --recode-INFO-all \
         --out ${OUTDIR}/minDP10_maxmiss0.8_filtInd
```

### ok now to compare some outputs! 

gstacks
```
Genotyped 178765 loci:
  effective per-sample coverage: mean=20.4x, stdev=11.1x, min=7.4x, max=47.6x
  mean number of sites per locus: 90.0
  a consistent phasing was found for 50503 of out 62535 (80.8%) diploid loci needing phasing
```

our prior run genotyped 178705 loci, so this is an improvement. the prior run also had a mean coverage of mean=20.3x. 

vcftools filtering (for the record i get this from the .sh.err file that gets output from stacks: so here, that command was `cat stacks_birostris_byc_own_pop_3034189.err' where the number code is the sbatch job #)

```
Parameters as interpreted:
	--vcf /projects/gatins/2025_Mobulid/birostris/byc_own_pop/populations.snps.vcf
	--recode-INFO-all
	--minDP 10
	--out /projects/gatins/2025_Mobulid/birostris/byc_own_pop/minDP10
	--recode

After filtering, kept 26 out of 26 Individuals
Outputting VCF file...
After filtering, kept 29075 out of a possible 29075 Sites
Run Time = 1.00 seconds

VCFtools - 0.1.16
(C) Adam Auton and Anthony Marcketta 2009

Parameters as interpreted:
	--vcf /projects/gatins/2025_Mobulid/birostris/byc_own_pop/minDP10.recode.vcf
	--recode-INFO-all
	--max-missing 0.8
	--out /projects/gatins/2025_Mobulid/birostris/byc_own_pop/minDP10_maxmiss0.8
	--recode

After filtering, kept 26 out of 26 Individuals
Outputting VCF file...
After filtering, kept 11826 out of a possible 29075 Sites
Run Time = 1.00 seconds

VCFtools - 0.1.16
(C) Adam Auton and Anthony Marcketta 2009

Parameters as interpreted:
	--vcf /projects/gatins/2025_Mobulid/birostris/byc_own_pop/minDP10_maxmiss0.8.recode.vcf
	--missing-indv
	--out /projects/gatins/2025_Mobulid/birostris/byc_own_pop/missingness

After filtering, kept 26 out of 26 Individuals
Outputting Individual Missingness
After filtering, kept 11826 out of a possible 11826 Sites
Run Time = 0.00 seconds

VCFtools - 0.1.16
(C) Adam Auton and Anthony Marcketta 2009

Parameters as interpreted:
	--vcf /projects/gatins/2025_Mobulid/birostris/byc_own_pop/minDP10_maxmiss0.8.recode.vcf
	--remove /projects/gatins/2025_Mobulid/birostris/byc_own_pop/remove_individuals.txt
	--recode-INFO-all
	--out /projects/gatins/2025_Mobulid/birostris/byc_own_pop/minDP10_maxmiss0.8_filtInd
	--recode

Excluding individuals in 'exclude' list
After filtering, kept 22 out of 26 Individuals
Outputting VCF file...
After filtering, kept 11826 out of a possible 11826 Sites
Run Time = 0.00 seconds
```

we lost a few loci here (down from 11843) but overall our results are looking really consistent. let's move on to the genetics analysis in R and make some plots! 

### transferring vcf to local R
I'm going to make a new folder on my desktop called birostris_byc_own_pop and store results there. if we end up going with this population map, I'll rename the folder to replace the current `/Users/madelineeppley/Desktop/manta/birostris` folder. 

``` bash
scp eppley.m@login.explorer.northeastern.edu:/projects/gatins/2025_Mobulid/birostris/byc_own_pop/minDP10_maxmiss0.8_filtInd.recode.vcf ~/Desktop/manta/birostris_byc_own_pop
```

### looking at some outputs
just taking a quick glance i'm not seeing any major differences between the first PCA and the PCA using results from grouping the bycatch sample with Peru. I'm going to retain these results just in case a reviewer asks for justificiation about why we grouped with Peru, but we should be able to stick with those same results from the initial run. 

### editing the DAPC
Here is the reviewer comment that I need to address: "Line 320 but relevant to other parts of the results as well: I would encourage the authors to be careful in implementation and interpretation of DAPC. There can be the appearance of structure, particularly in the cases in which populations are separated by very low FST separating populations. See recent (and semi-recent) published literature on the subject: 
https://onlinelibrary.wiley.com/doi/epdf/10.1111/1755-0998.13706
https://www.nature.com/articles/s41437-020-0348-2"

ok I went ahead and read through [Thia 2022](https://onlinelibrary.wiley.com/doi/epdf/10.1111/1755-0998.13706), which is the first paper listed by the reviewer to try to better understand DAPC analysis. DAPC is a hypothesis-driven test of a priori population groupings. The risk of DAPC is that overfitting can easily occur when too many explanatory PCs are retained, particularly in panmictic groups. 

In the previous version of the analysis, "Twenty principal components (PCs) and two discriminant axes (DAs) were retained for each
DAPC analysis." 

Here, let's say we are going to test a priori sampling locations as "populations" of m. birostris. 
- there were 3 true sampling locations (PERU, MEX, IND) and one bycatch location.
- accordingly, if we are going to include bycatch as it's own sampling location, the number of PCs retained to test this hypothesis for structure would be `n.pc = 4` and `n.da = 1 - n.pc`, which is 3. we would test this on the genetic dataset where we ran STACKS with the bycatch sample as it's own group.
- we could also test `n.pc = 3` and `n.da = 2` with the dataset where the bycatch sample was run through the STACKS program within the PERU group. this would account for 3 PCs retained for each of the three true sampling locations.

Here are the results: 

There are some very weak differences, and it does appear that our bycatch sample is separating just a bit from PERU. However, the main takeaway here is that the birostris popualtion is clearly panmictic across these sampling sites. It's worth noting that these findings do also align exactly with the [Humble et al 2023 manta paper](https://onlinelibrary.wiley.com/doi/epdf/10.1111/mec.17220), which investigated M. birostris and M. alfredi, finding population structure in alfredi but none in birostris. 
<img width="798" height="684" alt="image" src="https://github.com/user-attachments/assets/573e4e36-0c3d-492c-a3ad-a7f6bdbfae37" />


# Fst with the BYC as a separate pop 
actually shows more difference between IND and BYC sample - Fst increases from a max of 0.006 to 0.013
<img width="1593" height="1365" alt="image" src="https://github.com/user-attachments/assets/d4f3188a-300a-4342-af1a-0c25459ca783" />





