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
  --popmap /projects/gatins/2025_Mobulid/birostris \
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



