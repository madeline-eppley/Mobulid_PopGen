# Munkiana RADseq analysis

In this dataset, we have 29 RAD-sequenced samples from Peru (7), Mexico (12), and Nicaragua (10). 

Mobula munkiana, or the pygmy devil ray, is only found in coastal waters of the Eastern Pacific from Mexico to Peru. They form large schooling aggregates in the sea of Cortez and also breach! so cool :)

insane max vertical!!!

<img width="500" height="330" alt="image" src="https://github.com/user-attachments/assets/342a5286-49a7-4d35-8d0c-f845ae0a7b63" />

ok let's get into it 🫡

### move the pop map
move our munkiana pop map file into the new working dir with 
```bash
cp /projects/gatins/2025_Mobulid_UCSC/RAD_all_combined_bycatch/pop_map_munkiana /projects/gatins/2025_Mobulid/munkiana/pop_map_munkiana
```

Here's what we're working with! nice even distribution of samples across collection sites. 

```
mu_11	MEX
mu_12	MEX
mu_14	MEX
mu_19	MEX
mu_2	MEX
mu_21	MEX
mu_23	MEX
mu_25	MEX
mu_3	MEX
mu_4	MEX
mu_7	MEX
mu_8	MEX
NIC_2725_01	NIC
NIC_2726_01	NIC
NIC_2727_01	NIC
NIC_2728_01	NIC
NIC_2731_01	NIC
NIC_2731_02	NIC
NIC_2732_01	NIC
NIC_2736_01	NIC
NIC_2738_01	NIC
NIC_2740_01	NIC
NIC_6013_1	NIC
PER_DCLCW130_MM	PER
PER_DCLCW87_MM	PER
PER_DPPW019_MM	PER
PER_DPPW025_MM	PER
PER_DPPW033_MM	PER
PER_DPPW28_MM	PER
PER_DPPW98_MM	PER
```

## full pipeline
now let's execute our full stacks pipeline! we are using our standard settings here that we established after running the tarapacana dataset. however, i'm still going to go back and run n_opt scripts for plotting purposes. 

```bash
#!/bin/bash

#SBATCH --partition=lotterhos
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --time=48:00:00
#SBATCH --job-name=stacks_munkiana
#SBATCH --output=stacks_munkiana_%j.log
#SBATCH --error=stacks_munkiana_%j.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=eppley.m@northeastern.edu


export PATH=/projects/gatins/programs_explorer/stacks_2.68/bin:$PATH


echo "Sample list:"
cat /projects/gatins/2025_Mobulid/munkiana/pop_map_munkiana

echo "Starting denovo_map.pl"

denovo_map.pl \
  -m 3 \
  -M 3 \
  -n 2 \
  -T 32 \
  -o /projects/gatins/2025_Mobulid/munkiana \
  --popmap /projects/gatins/2025_Mobulid/munkiana/pop_map_munkiana \
  --samples /projects/gatins/2025_Mobulid_UCSC/RAD_all_combined_bycatch/trimmed90 \
  -X "populations:-r 0.8 -p 3 --min-maf 0.05 --write-single-snp --vcf --genepop --structure --fstats --hwe -t 30"


ls -l /projects/gatins/2025_Mobulid/munkiana

echo "populations module done"

echo "starting vcftools"
module load vcftools # this is just globally available on explorer, slay

# input VCF from stacks population run
INPUT_VCF="/projects/gatins/2025_Mobulid/munkiana/populations.snps.vcf"

OUTDIR="/projects/gatins/2025_Mobulid/munkiana"
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


### results 

gstacks output (cat gstacks.log)

```
Genotyped 213973 loci:
  effective per-sample coverage: mean=14.0x, stdev=7.1x, min=7.1x, max=30.8x
  mean number of sites per locus: 90.0
  a consistent phasing was found for 41043 of out 53231 (77.1%) diploid loci needing phasing
```

populations.log output (cat populations.log)

```
Removed 166929 loci that did not pass sample/population constraints from 213973 loci.
Kept 47044 loci, composed of 4235955 sites; 26189 of those sites were filtered, 6599 variant sites remained.
Mean genotyped sites per locus: 89.81bp (stderr 0.00).

Population summary statistics (more detail in populations.sumstats_summary.tsv):
  MEX: 11.268 samples per locus; pi: 0.29932; all/variant/polymorphic sites: 4225071/6599/6184; private alleles: 193
  NIC: 9.7148 samples per locus; pi: 0.28668; all/variant/polymorphic sites: 4225190/6599/5931; private alleles: 89
  PER: 6.8254 samples per locus; pi: 0.29818; all/variant/polymorphic sites: 4225021/6599/5611; private alleles: 23

Number of variable sites found to be significantly out of Hardy-Weinberg equilibrium (<0.05):
  MEX: 572
  NIC: 447
  PER: 318
Number of loci found to be significantly out of Hardy-Weinberg equilibrium (<0.05):
  MEX: 10773
  NIC: 10457
  PER: 10729
(more detail in populations.sumstats.tsv and populations.hapstats.tsv)

Population pair divergence statistics (more in populations.fst_summary.tsv and populations.phistats_summary.tsv):
  MEX-NIC: mean Fst: 0.031644; mean Phi_st: 0.010353; mean Fst': 0.0046198; mean Dxy: 0.003495
  MEX-PER: mean Fst: 0.032207; mean Phi_st: 0.0033874; mean Fst': 0.00098508; mean Dxy: 0.0035196
  NIC-PER: mean Fst: 0.034054; mean Phi_st: 0.0021441; mean Fst': -0.00074356; mean Dxy: 0.0035198
```

we did end up removing three individuals from the Nicaragua population. I will need to check again what the pop sizes are for each sampling site. 

```
(base) [eppley.m@explorer-01 munkiana]$ cat remove_individuals.txt 
INDV
NIC_2726_01
NIC_2731_01
NIC_2731_02
```


depths of coverage for this species are on average lower than some of the other species, but still should be fine: 
```
Depths of Coverage for Processed Samples:
mu_11; loci assembled: 127301; depth: 19.65x; max: 958x; number of reads: 2487282 (81.4%)
mu_12; loci assembled: 127064; depth: 18.65x; max: 1273x; number of reads: 2356285 (78.2%)
mu_14; loci assembled: 91217; depth: 7.64x; max: 485x; number of reads: 693002 (69.8%)
mu_19; loci assembled: 101187; depth: 7.82x; max: 436x; number of reads: 787444 (74.7%)
mu_2; loci assembled: 106319; depth: 7.89x; max: 317x; number of reads: 834841 (73.0%)
mu_21; loci assembled: 123282; depth: 22.79x; max: 1077x; number of reads: 2792707 (81.5%)
mu_23; loci assembled: 91954; depth: 7.51x; max: 447x; number of reads: 687252 (74.6%)
mu_25; loci assembled: 107813; depth: 8.74x; max: 573x; number of reads: 938023 (81.1%)
mu_3; loci assembled: 128901; depth: 29.48x; max: 1152x; number of reads: 3777805 (80.6%)
mu_4; loci assembled: 127861; depth: 22.32x; max: 1030x; number of reads: 2836042 (81.4%)
mu_7; loci assembled: 127106; depth: 19.39x; max: 1142x; number of reads: 2451037 (81.0%)
mu_8; loci assembled: 96419; depth: 7.64x; max: 389x; number of reads: 733463 (70.9%)
NIC_2725_01; loci assembled: 124082; depth: 12.06x; max: 544x; number of reads: 1489906 (82.7%)
NIC_2726_01; loci assembled: 76712; depth: 7.06x; max: 292x; number of reads: 539145 (77.8%)
NIC_2727_01; loci assembled: 68233; depth: 7.67x; max: 281x; number of reads: 521264 (78.8%)
NIC_2728_01; loci assembled: 139444; depth: 27.19x; max: 1024x; number of reads: 3772052 (86.0%)
NIC_2731_01; loci assembled: 63703; depth: 7.00x; max: 307x; number of reads: 443602 (78.0%)
NIC_2731_02; loci assembled: 60125; depth: 7.15x; max: 179x; number of reads: 428304 (79.1%)
NIC_2732_01; loci assembled: 104168; depth: 12.31x; max: 553x; number of reads: 1275975 (83.6%)
NIC_2736_01; loci assembled: 103266; depth: 10.83x; max: 450x; number of reads: 1112540 (82.6%)
NIC_2738_01; loci assembled: 123791; depth: 12.33x; max: 507x; number of reads: 1518993 (84.8%)
NIC_2740_01; loci assembled: 134957; depth: 24.19x; max: 682x; number of reads: 3246446 (85.7%)
NIC_6013_1; loci assembled: 64973; depth: 8.01x; max: 407x; number of reads: 518368 (79.4%)
PER_DCLCW130_MM; loci assembled: 117700; depth: 13.54x; max: 399x; number of reads: 1584674 (85.2%)
PER_DCLCW87_MM; loci assembled: 120286; depth: 16.01x; max: 495x; number of reads: 1917031 (87.0%)
PER_DPPW019_MM; loci assembled: 118409; depth: 15.63x; max: 474x; number of reads: 1840668 (84.7%)
PER_DPPW025_MM; loci assembled: 101625; depth: 8.14x; max: 318x; number of reads: 823244 (80.4%)
PER_DPPW033_MM; loci assembled: 106851; depth: 9.04x; max: 249x; number of reads: 961526 (83.4%)
PER_DPPW28_MM; loci assembled: 110690; depth: 10.26x; max: 293x; number of reads: 1130484 (85.0%)
PER_DPPW98_MM; loci assembled: 125316; depth: 12.97x; max: 431x; number of reads: 1613917 (79.4%)
```


