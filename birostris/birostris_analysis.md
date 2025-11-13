# Birostris RADseq analysis

In this dataset, we have 36 RAD-sequenced samples from India, Peru, Mexico, and one bycatch sample from the high seas.

Mobula birostris, or the giant oceanic manta ray, is the largest ray in the world :o 

It is an endangered species in tropical and temperate waters across the world. let's get into some genetics!!!

<img width="554" height="554" alt="image" src="https://github.com/user-attachments/assets/ec6d715f-909b-4505-8c5b-4aa676f1c872" />


### fix the pop map
first things first let's fix this pop map: 
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


## full pipeline
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
  -o /projects/gatins/2025_Mobulid/birostris/ \
  --popmap /projects/gatins/2025_Mobulid/birostris \
  --samples /projects/gatins/2025_Mobulid_UCSC/RAD_all_combined_bycatch/trimmed90 \
  -X "populations:-r 0.8 -p 1 --min-maf 0.05 --write-single-snp --vcf --genepop --structure --fstats --hwe -t 30"


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
