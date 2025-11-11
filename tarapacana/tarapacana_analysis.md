# Tarapacana RADseq analysis

we are working with bycatch samples from the high seas for the species Mobula tarapacana. We have 15 RAD-sequenced samples that were collected between ~2018-2022. 

first things first! here's the [critically endangered](https://www.iucnredlist.org/es/species/60199/279077133) Mobula tarapacana, or the sicklefin devil ray. wide-ranging tropical/temperate species.

<img width="500" height="412" alt="image" src="https://github.com/user-attachments/assets/93c78dfe-9f71-4f1d-b10d-b7eeeed9556a" />

### ok now for the fun stuff!

#### Stacks denovo map
files uploaded from melissa and located here:
```bash
/projects/gatins/2025_Mobulid_UCSC/RAD_all_combine_bycatch
```

then we have subfolders for each species, including tarapacana

so this species is unique because we have all high-seas samples. we have lat/long info for 14/15 samples (BYC_RMT_59 is missing)

`pop_map_tarapacana`

```bash
BYC_RMB_56	BYC
BYC_RMB_57	BYC
BYC_RMM_30	BYC
BYC_RMO_45	BYC
BYC_RMT_04	BYC
BYC_RMT_06	BYC
BYC_RMT_07	BYC
BYC_RMT_27	BYC
BYC_RMT_28	BYC
BYC_RMT_29	BYC
BYC_RMT_46	BYC
BYC_RMT_49	BYC
BYC_RMT_59	BYC
BYCI_RMT_69	BYC
BYCI_RMT_71	BYC
```

Let's look at some of the lat/long info from the excel sheet:


| Sample | Latitude | Longitude |
|-----------|----------|-----------|
| BYC_RMM_30 | 2°43'00.0"N | 91°59'00.0"W |
| BYC_RMB_57 | 1°26'00.0"N | 82°26'00.0"W |
| BYC-RMB-56 | 14°03'00.0"N | 82°34'00.0"W |
| BYC_RMO_45 | 2°13'00.0"N | 85°19'00.0"W |
| BYC_RMT_04 | 3°49'12.0"N | 110°04'48.0"W |
| BYC_RMT_06 | 4°16'59.9"N | 110°04'59.9"W |
| BYC_RMT_07 | 0°54'00.0"S | 92°13'48.0"W |
| BYC_RMT_27 | 2°39'00.0"N | 84°37'00.1"W |
| BYC_RMT_69 | 4°19'00.1"N | 144°24'00.0"W |
| BYC_RMT_71 | 4°36'00.0"S | 179°52'48.0"E |
| BYC_RMT_49 | 2°24'00.0"N | 90°04'48.0"W |
| BYC_RMT_28 | 3°21'00.0"N | 92°06'00.0"W |
| BYC_RMT_29 | 0°21'00.0"N | 82°19'58.8"W |
| BYC_RMT_46 | 1°46'12.0"N | 84°48'00.0"W |


<img width="1447" height="572" alt="Screenshot 2025-11-10 at 10 09 19 AM" src="https://github.com/user-attachments/assets/44fcd6ce-b5f7-419b-9e82-4acdfe730970" />


based on this output, I think a good starting place would be 3 for populations, so I made a new pop map file

this file has an eastern pacific (EST), central pacific (CTR), and western pacific group structure (WST). we are missing for BYC_RMT_59, but I assigned it eastern to start and we will go from there. 

```bash
cat > tarapacana_pop_map_3 << 'EOL'
BYC_RMB_56 EST
BYC_RMB_57 EST
BYC_RMM_30 EST
BYC_RMO_45 EST
BYC_RMT_04 CTR
BYC_RMT_06 CTR
BYC_RMT_07 EST
BYC_RMT_27 EST
BYC_RMT_28 EST
BYC_RMT_29 EST
BYC_RMT_46 EST
BYC_RMT_49 EST
BYC_RMT_59 EST
BYCI_RMT_69 CTR
BYCI_RMT_71 WST
EOL
```

I want to get an initial version of the denovo map up and working, and then we can go from there with parameter testing. 

Inital test run of this denovo map script works! 

```bash
#!/bin/bash

#SBATCH --partition=lotterhos
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --time=48:00:00
#SBATCH --job-name=stacks_tarapacana
#SBATCH --output=stacks_tarapacana_%j.log
#SBATCH --error=stacks_tarapacana_%j.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=eppley.m@northeastern.edu


export PATH=/projects/gatins/programs_explorer/stacks_2.68/bin:$PATH



SAMPLE_COUNT=$(wc -l < /projects/gatins/2025_Mobulid_UCSC/RAD_all_combined_bycatch/pop_map_tarapacana)
echo "Total num of samples: $SAMPLE_COUNT"

echo "Sample list:"
cat /projects/gatins/2025_Mobulid_UCSC/RAD_all_combined_bycatch/pop_map_tarapacana

echo "Starting denovo_map.pl"

denovo_map.pl \
  -m 3 \
  -M 3 \
  -n 4 \
  -T 32 \
  -o /projects/gatins/2025_Mobulid/tarapacana \
  --popmap /projects/gatins/2025_Mobulid_UCSC/RAD_all_combined_bycatch/pop_map_tarapacana \
  --samples /projects/gatins/2025_Mobulid_UCSC/RAD_all_combined_bycatch/trimmed90



ls -l /projects/gatins/2025_Mobulid/tarapacana
```

## optimizing parameters
ok time to get serious. we've got stacks up and running and now we should focus on the parameter optimization. We can determine the optimal params from the highest yield loci at the end of the pipeline. 

#### stacks params
**Options:**
`-m` — Minimum depth of coverage required to create a stack (default 3).  
`-M` — number of mismatches allowed between stacks within individuals (for ustacks).  
`-n` — number of mismatches allowed between stacks between individuals (for cstacks).  
`-samples [path]` — specify a path to the directory of samples (samples will be read from population map).  
`--popmap [path]` — path to a population map file (format is "[name] TAB [pop]", one sample per line).  
`-o [path]` — path to write pipeline output files.  
`-p`,`--min-populations` — minimum number of populations a locus must be present in to process a locus (for populations; default: 1). 
`-r`,`--min-samples-per-pop` — minimum percentage of individuals in a population required to process a locus for that population (for populations; default: 0). 

#### Remy and I agree that -m 3 and -M 3 are default, let's keep those and instead vary parameter -n
for some further proof beyond the BSB dataset, here's what [Paris et al 2017](https://besjournals.onlinelibrary.wiley.com/doi/epdf/10.1111/2041-210X.12775) says about stacks parameter space across 4 different datasets: "We recommend setting m high enough only to deny errantreads the status of a putative allele. We therefore do not sug-gest using a value for m > 5 and we have demonstrated herethat the default value of m3 was favourable for all test datasets"

tldr: -m 3 is default and should be good here

```bash
(base) [eppley.m@explorer-02 tarapacana]$ cat tarapacana_opt_n.sh 
#!/bin/bash

#SBATCH --partition=lotterhos
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --time=48:00:00
#SBATCH --job-name=stacks_tarapacana
#SBATCH --output=stacks_tarapacana_%j.log
#SBATCH --error=stacks_tarapacana_%j.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=eppley.m@northeastern.edu

export PATH=/projects/gatins/programs_explorer/stacks_2.68/bin:$PATH

# new dir for our opt folder
mkdir -p /projects/gatins/2025_Mobulid/tarapacana/opt

# parameters to test for n
params="
2
3
4
5
6
7"

# now denovo_map.pl for each n value
for p in $params
do
    echo "Running with n = ${p}"
    mkdir -p /projects/gatins/2025_Mobulid/tarapacana/opt/n${p}
    
    denovo_map.pl \
        -m 3 \
        -M 3 \
        -n ${p} \
        -T 32 \
        -o /projects/gatins/2025_Mobulid/tarapacana/opt/n${p} \
        --popmap /projects/gatins/2025_Mobulid_UCSC/RAD_all_combined_bycatch/pop_map_tarapacana \
        --samples /projects/gatins/2025_Mobulid_UCSC/RAD_all_combined_bycatch/trimmed90
done
```

### summary table 

using `cat populations.hapstats.tsv | grep -v "^#" | cut -f 1 | uniq | wc -l` in each folder of n2, n3, n4, etc.

n2 = 81560 loci
n3 = 81363 loci
n4 = 81047 loci
n5 = 80505 loci


## testing values of p (p = pops)
now we want to get into population testing since we don't have pre-existing expectations about # of pops in the data set. 
for our optimal catalog with `-m 3 -M 3 -n opt`, we are going to cycle through the populations with different pop maps that each have different groups of samples based on collection location. 

#### let's organize 🫡
step 1: move a copy of the p=1 og pop map into the new pop map folder and rename it to be pop_map_1 to make my life easier

`cp /projects/gatins/2025_Mobulid_UCSC/RAD_all_combined_bycatch/pop_map_tarapacana /projects/gatins/2025_Mobulid/tarapacana/pop_maps/tarapacana_pop_map_1`

and now repeat to do the same with our new n=3 pop map that we made earlier based on the high seas sampling locations

`mv /projects/gatins/2025_Mobulid/tarapacana/tarapacana_pop_map_3 /projects/gatins/2025_Mobulid/tarapacana/pop_maps/`

let's make a pop map file for p=2 by combining central with western:

```bash
cat > /projects/gatins/2025_Mobulid/tarapacana/pop_maps/tarapacana_pop_map_2 << 'EOF'
BYC_RMB_56 EST
BYC_RMB_57 EST
BYC_RMM_30 EST
BYC_RMO_45 EST
BYC_RMT_04 WST
BYC_RMT_06 WST
BYC_RMT_07 EST
BYC_RMT_27 EST
BYC_RMT_28 EST
BYC_RMT_29 EST
BYC_RMT_46 EST
BYC_RMT_49 EST
BYC_RMT_59 EST
BYCI_RMT_69 WST
BYCI_RMT_71 WST
EOF
```

so at this point, p=2 and p=3 are kind of the end of my educated "guessing" based on geography. we will also examine the PCA with p=1 to check for signs of population structure. the last thing we can do is make every individual their OWN pop, and check the PCA there too. so let's make that pop map where every ind is its own pop:

```bash
cat > /projects/gatins/2025_Mobulid/tarapacana/pop_maps/pop_map_15 << 'EOF'
BYCI_RMT_69 1
BYCI_RMT_71 2
BYC_RMB_56 3
BYC_RMB_57 4
BYC_RMM_30 5
BYC_RMO_45 6
BYC_RMT_04 7
BYC_RMT_06 8
BYC_RMT_07 9
BYC_RMT_27 10
BYC_RMT_28 11
BYC_RMT_29 12
BYC_RMT_46 13
BYC_RMT_49 14
BYC_RMT_59 15
EOF
```


### opt params populations test
the n opt is 2 based on the loci yield. let's get the populations script running with this n value for our p=1 map. 

#### gstacks output

```bash
Genotyped 174775 loci:
  effective per-sample coverage: mean=26.4x, stdev=10.5x, min=10.2x, max=54.5x
  mean number of sites per locus: 90.1
  a consistent phasing was found for 73361 of out 81544 (90.0%) diploid loci needing phasing
```


#### cstacks output

```bash
BEGIN cov_per_sample
# Depths of Coverage for Processed Samples
sample	loci assembled	depth of cov	max cov	number reads incorporated	% reads incorporated
BYC_RMB_56	133278	19.61	335	2594594	85.3
BYC_RMB_57	132472	26.76	376	3516746	85.7
BYC_RMM_30	131071	20.83	491	2709256	85.7
BYC_RMO_45	133560	29.78	673	3944421	85.7
BYC_RMT_04	131613	23.26	480	3036583	86.7
BYC_RMT_06	129428	19.52	503	2507191	86.6
BYC_RMT_07	135653	52.05	617	6981223	85.9
BYC_RMT_27	131345	21.80	417	2839859	86.2
BYC_RMT_28	131301	27.91	428	3633508	86.7
BYC_RMT_29	123942	11.87	263	1459511	83.9
BYC_RMT_46	134070	33.41	566	4440150	85.8
BYC_RMT_49	133206	29.52	614	3898435	85.7
BYC_RMT_59	123544	10.05	308	1232876	80.1
BYCI_RMT_69	133092	37.07	702	4889100	87.0
BYCI_RMT_71	133038	20.75	451	2738769	85.6
END cov_per_sample
```

## running populations

```bash
#!/bin/bash

#SBATCH --partition=lotterhos
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --time=48:00:00
#SBATCH --job-name=tarapacana_pops
#SBATCH --output=tarapacana_pops_%j.log
#SBATCH --error=tarapacana_pops_%j.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=eppley.m@northeastern.edu

export PATH=/projects/gatins/programs_explorer/stacks_2.68/bin:$PATH


mkdir -p /projects/gatins/2025_Mobulid/tarapacana/populations/n2_p1

populations \
  -P /projects/gatins/2025_Mobulid/tarapacana/opt/n2 \
  -M /projects/gatins/2025_Mobulid/tarapacana/pop_maps/tarapacana_pop_map_1 \
  -r 0.80 \
  -p 1 \
  --min-maf 0.05 \
  --write-single-snp \
  --vcf \
  --genepop \
  --structure \
  --fstats \
  --hwe \
  -O /projects/gatins/2025_Mobulid/tarapacana/populations/n2_p1 \
  -t 30

```
### MAF and LD filtered loci
after running populations I counted loci with:

```bash
grep -v "^#" /projects/gatins/2025_Mobulid/tarapacana/populations/n2_p1/populations.hapstats.tsv | wc -l
```
#### 73660 loci remain



## running vcftools for filtering loci
we need to make a rigorous plan for how we will deal with linkage, minor alleles, and missingness. all of these can bias population structure results. 

first off, `--write_single_snp` that we just did filtered for linkage. from the stacks manual - " Since STRUCTURE does not want linked loci, you will typically also supply the --write_single_snp flag so that you do not get more than one SNP per RAD locus (SNPs at the same RAD locus are almost certainly linked)." 

let's use vcftools to filter for loci that don't meet a minimum depth, loci that aren't present in most individuals, and individuals with overall poor coverage. 


### filtering results

##### individuals
unfortunately it looks like our BYC_RMT_59 sample didn't pass our filtering methods. we'll have to remove it. 

```bash
(base) [eppley.m@explorer-02 vcftools_filtered]$ cd n2_p1
(base) [eppley.m@explorer-02 n2_p1]$ ls
minDP10_maxmiss0.8_filtInd.recode.vcf  minDP10_maxmiss0.8.recode.vcf  minDP10.recode.vcf  missingness.imiss  remove_individuals.txt
(base) [eppley.m@explorer-02 n2_p1]$ cat remove_individuals.txt 
INDV
BYC_RMT_59
```


```bash
#!/bin/bash

#SBATCH --partition=lotterhos
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --job-name=vcftools
#SBATCH --output=vcftools_filter_%j.log
#SBATCH --error=vcftools_filter_%j.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=eppley.m@northeastern.edu

module load vcftools # this is just globally available on explorer, slay

# input VCF from stacks population run
INPUT_VCF="/projects/gatins/2025_Mobulid/tarapacana/populations/n2_p1/populations.snps.vcf"

OUTDIR="/projects/gatins/2025_Mobulid/tarapacana/vcftools_filtered/n2_p1"
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

### getting vcfs into R locally
transfer files with scp and the following path:
`scp eppley.m@login.explorer.northeastern.edu:/projects/gatins/2025_Mobulid/tarapacana/vcftools_filtered/n2_p1/minDP10_maxmiss0.8_filtInd.recode.vcf ~/Users/madelineeppley/Desktop/manta`

## pop structure time!!!
let's start with a PCA

```R
library(vcfR)
library(adegenet)

vcf <- read.vcfR("/Users/madelineeppley/Desktop/manta/minDP10_maxmiss0.8_filtInd.recode.vcf")
# we have 48,350 SNPs remaining

# convert to genlight
gl <- vcfR2genlight(vcf)

pca <- glPca(gl, nf=3)
#scatter(pca) # don't like these labels

pca_df <- data.frame(
  PC1 = pca$scores[,1],
  PC2 = pca$scores[,2],
  Sample = indNames(gl))

# plot again with better labels
ggplot(pca_df, aes(x=PC1, y=PC2)) +
  geom_point(size=2, color="blue") +
  geom_text_repel(aes(label=Sample),
                  size=2.5,
                  max.overlaps=50,
                  segment.size=0.3) +
  theme_minimal() +
  labs(title="tarapacana n2_p1",
       x=paste0("PC1 (", round(pca$eig[1]/sum(pca$eig)*100,1), "%)"),
       y=paste0("PC2 (", round(pca$eig[2]/sum(pca$eig)*100,1), "%)"))
```

<img width="1671" height="1224" alt="image" src="https://github.com/user-attachments/assets/1e180173-5c98-4614-84a0-4b42f79a7cd0" />


weird happenings with BYC_RMB_56 and BYC_RMB_57!

let's go back to what we know about them. they were collected very close together off the coast of ecuador, but in a similar area to some other samples. they were also collected just one day apart ... however one is male and one is female, so we shouldn't have a sex-related effect

<img width="1436" height="659" alt="Screenshot 2025-11-11 at 8 13 09 AM" src="https://github.com/user-attachments/assets/48e7fb18-4455-4dcd-91f3-ef99f5e740d5" />

collection date, species, sex, lat, long, sample: 

7/14/21	RMB	H	0143N	08234W	BYC_RMB_56

7/15/21	RMB	M	0126N	08226W	BYC_RMB_57


possible options: 
- we have relatedness between these two samples (BYC_RMB_56 and BYC_RMB_57) such as offspring/parent, sibling pair
- these are a different species? i.e. cryptic speciation
- these samples have an artifical effect like high missingness 
- we have bad SNP filtering and should go back to the drawing board.

#### missingness at the individual level
Our BYC_RMB_56 and BYC_RMB_57 individuals have low missingness ~4% for 56 and ~.07% for 57, so it's not that! However BYC_RMT_29 is borderline ...

```R
# check missingness per ind
gt <- extract.gt(vcf, element = "GT")
ind_missing <- apply(gt, 2, function(x) mean(is.na(x) | x == "./."))
snp_missing <- apply(gt, 1, function(x) mean(is.na(x) | x == "./."))
summary(ind_missing)
summary(snp_missing)

sort(ind_missing, decreasing = TRUE)[1:14]
```

output
```
> sort(ind_missing, decreasing = TRUE)[1:14]
 BYC_RMT_29  BYC_RMT_06  BYC_RMB_56 BYCI_RMT_71  BYC_RMM_30  BYC_RMT_27 
0.338572906 0.077642192 0.040972079 0.031023785 0.030568769 0.030361944 
 BYC_RMT_04  BYC_RMT_28  BYC_RMO_45  BYC_RMB_57  BYC_RMT_46  BYC_RMT_07 
0.022647363 0.016504654 0.007797311 0.007487073 0.006866598 0.006287487 
 BYC_RMT_49 BYCI_RMT_69 
0.006204757 0.005729059 
```

overall the dataset looks pretty fire
<img width="1671" height="1224" alt="image" src="https://github.com/user-attachments/assets/bfb12201-5987-408b-817f-ac05b228b1a1" />
<img width="1671" height="1224" alt="image" src="https://github.com/user-attachments/assets/3194e64a-1684-4d86-a68c-c2f9b2c43d41" />

#### cryptic speciation ...??
I will note that these two samples were identified as species RMB (birostris) by the fishers who collected which I find odd. however the STRUCTURE results from the misidentification report clearly group those 2 samples with RMT with very high support (100 bootstrap values on phylogenetic tree). It seems possible to me that the scale of differentiation between the species in the phylogenetics study could obstruct some finer scale cryptic speciation within species?

<img width="468" height="508" alt="image" src="https://github.com/user-attachments/assets/6e8cf2b2-876d-451f-8e95-7dbea792f300" />

### relatedness 🙂‍↕️
we have a high level of relatedness! 0.5 = an exact clone or identical twins. Here we have a value of 0.498 between RMB_56 and RMB_57. I think it's highly likely this was the same individual sampled twice. Manta rays usually have just 1 offspring at a time, and the chances of monozygotic twins are super low. HOWEVER the one super weird thing here is that the sampling sheet says 56 is a Female and 57 is a Male. 

<img width="2013" height="1596" alt="image" src="https://github.com/user-attachments/assets/20e7cccd-55a7-4187-b24e-39e65efafc27" />

```R
IBD analysis (KING method of moment) on genotypes:
Excluding 319 SNPs (monomorphic: TRUE, MAF: NaN, missing rate: NaN)
    # of samples: 14
    # of SNPs: 48,031
    using 1 thread
No family is specified, and all individuals are treated as singletons.
Relationship inference in the presence of population stratification.
KING IBD:    the sum of all selected genotypes (0,1,2) = 1034229
```

```R
head(rel_df, 10)
    ind1 ind2     kinship
2      B    A 0.498100760
82     L    F 0.008920038
34     F    C 0.006222788
38     J    C 0.005227273
52     J    D 0.005096953
137    K    J 0.005090748
80     J    F 0.004961348
66     J    E 0.002706775
94     J    G 0.002628375
138    L    J 0.002402253
```

### snmf ancestry coefficients
with both samples in, testing K between 1:5, our lowest cross-entropy value is K=2

```R
$crossEntropy
         K = 1     K = 2    K = 3    K = 4    K = 5
min  0.8863364 0.8454335 1.056315 1.248850 1.401853
mean 0.8906090 0.8528392 1.075903 1.271387 1.430143
max  0.8975796 0.8630781 1.099681 1.294355 1.464540
```

<img width="1938" height="1416" alt="image" src="https://github.com/user-attachments/assets/b4b31bf3-fa82-4047-b358-4756cd164664" />


## back to the filtering & doing some sample elimination
changes we can make this time: use a more strict MAF 0.01; run the whole denovo and populations pipeline at the same time; remove 1 sample between 56 and 57; remove BOTH samples 56 and 57. 


### remove 1 sample
started a new directory called /elim/rmv_1 and removed sample 56, which had higher missingness than sample 57 (~4% vs 0.7%) 

`cp /projects/gatins/2025_Mobulid_UCSC/RAD_all_combined_bycatch/pop_map_tarapacana \
   /projects/gatins/2025_Mobulid/tarapacana/elim/rmv_1/`


```bash
(base) [eppley.m@explorer-02 rmv_1]$ cat pop_map_rmv1_tarapacana 
BYC_RMB_57	BYC
BYC_RMM_30	BYC
BYC_RMO_45	BYC
BYC_RMT_04	BYC
BYC_RMT_06	BYC
BYC_RMT_07	BYC
BYC_RMT_27	BYC
BYC_RMT_28	BYC
BYC_RMT_29	BYC
BYC_RMT_46	BYC
BYC_RMT_49	BYC
BYC_RMT_59	BYC
BYCI_RMT_69	BYC
BYCI_RMT_71	BYC
```

## full pipeline!

```bash
(base) [eppley.m@explorer-02 tarapacana]$ cat one_removed_tarapacana_pipeline.sh 
#!/bin/bash

#SBATCH --partition=lotterhos
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --time=48:00:00
#SBATCH --job-name=stacks_tarapacana
#SBATCH --output=elim1_stacks_tarapacana_%j.log
#SBATCH --error=elim1_stacks_tarapacana_%j.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=eppley.m@northeastern.edu


export PATH=/projects/gatins/programs_explorer/stacks_2.68/bin:$PATH


echo "Sample list:"
cat /projects/gatins/2025_Mobulid/tarapacana/elim/rmv_1/pop_map_rmv1_tarapacana

echo "Starting denovo_map.pl"

denovo_map.pl \
  -m 3 \
  -M 3 \
  -n 2 \
  -T 32 \
  -o /projects/gatins/2025_Mobulid/tarapacana/elim/rmv_1 \
  --popmap /projects/gatins/2025_Mobulid/tarapacana/elim/rmv_1/pop_map_rmv1_tarapacana \
  --samples /projects/gatins/2025_Mobulid_UCSC/RAD_all_combined_bycatch/trimmed90 \
  -X "populations:-r 0.8 -p 1 --min-maf 0.05 --write-single-snp --vcf --genepop --structure --fstats --hwe -t 30"


ls -l /projects/gatins/2025_Mobulid/tarapacana/elim/rmv_1

echo "populations module done"

echo "starting vcftools"
module load vcftools # this is just globally available on explorer, slay

# input VCF from stacks population run
INPUT_VCF="/projects/gatins/2025_Mobulid/tarapacana/elim/rmv_1/populations.snps.vcf"

OUTDIR="/projects/gatins/2025_Mobulid/tarapacana/elim/rmv_1"
mkdir -p ${OUTDIR}

# filter by minimum depth per genotype (minDP = 10)
vcftools --vcf ${INPUT_VCF} \
         --minDP 10 \
         --recode --recode-INFO-all \
         --out ${OUTDIR}/minDP10_elim1

# filter for sites present in >= 80% of individuals
vcftools --vcf ${OUTDIR}/minDP10_elim1.recode.vcf \
         --max-missing 0.8 \
         --recode --recode-INFO-all \
         --out ${OUTDIR}/minDP10_maxmiss0.8_elim1

# remove individuals with >40% missing data
# 1: compute missingness per individual
vcftools --vcf ${OUTDIR}/minDP10_maxmiss0.8_elim1.recode.vcf \
         --missing-indv \
         --out ${OUTDIR}/missingness_elim1

# 2: generate a list of individuals to remove
awk '$5 > 0.4 {print $1}' ${OUTDIR}/missingness_elim1.imiss > ${OUTDIR}/remove_individuals_elim1.txt

# now filter out individuals
vcftools --vcf ${OUTDIR}/minDP10_maxmiss0.8_elim1.recode.vcf \
         --remove ${OUTDIR}/remove_individuals_elim1.txt \
         --recode --recode-INFO-all \
         --out ${OUTDIR}/minDP10_maxmiss0.8_filtInd_elim1
```

and now repeat for removing BOTH samples in our test

```bash
#!/bin/bash

#SBATCH --partition=lotterhos
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --time=48:00:00
#SBATCH --job-name=stacks_tarapacana
#SBATCH --output=elim2_stacks_tarapacana_%j.log
#SBATCH --error=elim2_stacks_tarapacana_%j.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=eppley.m@northeastern.edu


export PATH=/projects/gatins/programs_explorer/stacks_2.68/bin:$PATH


echo "Sample list:"
cat /projects/gatins/2025_Mobulid/tarapacana/elim/rmv_2/pop_map_rmv2_tarapacana

echo "Starting denovo_map.pl"

denovo_map.pl \
  -m 3 \
  -M 3 \
  -n 2 \
  -T 32 \
  -o /projects/gatins/2025_Mobulid/tarapacana/elim/rmv_2 \
  --popmap /projects/gatins/2025_Mobulid/tarapacana/elim/rmv_2/pop_map_rmv2_tarapacana \
  --samples /projects/gatins/2025_Mobulid_UCSC/RAD_all_combined_bycatch/trimmed90 \
  -X "populations:-r 0.8 -p 1 --min-maf 0.05 --write-single-snp --vcf --genepop --structure --fstats --hwe -t 30"


ls -l /projects/gatins/2025_Mobulid/tarapacana/elim/rmv_2

echo "populations module done"

echo "starting vcftools"
module load vcftools # this is just globally available on explorer, slay

# input VCF from stacks population run
INPUT_VCF="/projects/gatins/2025_Mobulid/tarapacana/elim/rmv_2/populations.snps.vcf"

OUTDIR="/projects/gatins/2025_Mobulid/tarapacana/elim/rmv_2"
mkdir -p ${OUTDIR}

# filter by minimum depth per genotype (minDP = 10)
vcftools --vcf ${INPUT_VCF} \
         --minDP 10 \
         --recode --recode-INFO-all \
         --out ${OUTDIR}/minDP10_elim2

# filter for sites present in >= 80% of individuals
vcftools --vcf ${OUTDIR}/minDP10_elim2.recode.vcf \
         --max-missing 0.8 \
         --recode --recode-INFO-all \
         --out ${OUTDIR}/minDP10_maxmiss0.8_elim2

# remove individuals with >40% missing data
# 1: compute missingness per individual
vcftools --vcf ${OUTDIR}/minDP10_maxmiss0.8_elim2.recode.vcf \
         --missing-indv \
         --out ${OUTDIR}/missingness_elim2

# 2: generate a list of individuals to remove
awk '$5 > 0.4 {print $1}' ${OUTDIR}/missingness_elim2.imiss > ${OUTDIR}/remove_individuals_elim2.txt

# now filter out individuals
vcftools --vcf ${OUTDIR}/minDP10_maxmiss0.8_elim2.recode.vcf \
         --remove ${OUTDIR}/remove_individuals_elim2.txt \
         --recode --recode-INFO-all \
         --out ${OUTDIR}/minDP10_maxmiss0.8_filtInd_elim2
```

## eliminate 1 of the highly related samples results

The total number of samples: 13 
The total number of SNPs: 43369 
244 SNPs - monomorphic: TRUE

PCA - nice, not seeing a relatedness effect here 

<img width="1569" height="1416" alt="image" src="https://github.com/user-attachments/assets/62d71622-05bf-40e6-8491-592bf60a7de9" />

missingness - nice results here, really complete data

```R
> summary(ind_missing_elim1)
    Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
0.003159 0.004289 0.008186 0.034521 0.018884 0.299361 
> summary(snp_missing_elim1)
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
0.00000 0.00000 0.00000 0.03452 0.07692 0.15385 
> 
> sort(ind_missing_elim1, decreasing = TRUE)[1:13]
 BYC_RMT_29  BYC_RMT_06  BYC_RMM_30 BYCI_RMT_71  BYC_RMT_27  BYC_RMT_04 
0.299361295 0.047683829 0.019668427 0.018884457 0.017132053 0.012889391 
 BYC_RMT_28  BYC_RMB_57  BYC_RMO_45  BYC_RMT_46  BYC_RMT_07  BYC_RMT_49 
0.008185570 0.005280269 0.004703821 0.004288778 0.004035140 0.003504808 
BYCI_RMT_69 
0.003158938
```

relatedness slay

<img width="1569" height="1416" alt="image" src="https://github.com/user-attachments/assets/d60c190b-85f4-4354-a15e-0c2064477cf2" />
