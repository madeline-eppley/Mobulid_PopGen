# Munkiana RADseq analysis

In this dataset, we have 29 RAD-sequenced samples from Peru (7), Mexico (12), and Nicaragua (10). 

Mobula munkiana, or the pygmy devil ray, is only found in coastal waters of the Eastern Pacific from Mexico to Peru. They form large schooling aggregates in the sea of Cortez and also breach! so cool :)

insane max vertical!!!
<img width="2048" height="1362" alt="image" src="https://github.com/user-attachments/assets/342a5286-49a7-4d35-8d0c-f845ae0a7b63" />

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
