## this is how to filter for long-range LD using PLINK
caveat: I realized after running this that the `write-single-snp` output from STACKS functionally saves 1 locus per chromosome. As a result, PLINK does not work becuase it is using a sliding window to filter WITHIN chromosomes. 
However, I will still upload my code here in case it is useful in the future. I also converted VCF -> PLINK, which could be helpful for other applications. 

the output from the PLINK run
```bash
(base) [eppley.m@explorer-02 birostris]$ cat plink_birostris_2951554.err

VCFtools - 0.1.16
(C) Adam Auton and Anthony Marcketta 2009

Parameters as interpreted:
	--vcf minDP10_maxmiss0.8_filtInd.recode.vcf
	--out birostris_minDP10_maxmiss0.8_filtInd
	--plink

After filtering, kept 22 out of 22 Individuals
Writing PLINK PED and MAP files ... 
Done.
After filtering, kept 11843 out of a possible 11843 Sites
```

```bash
(base) [eppley.m@explorer-02 birostris]$ cat plink2filter_birostris.sh 
#!/bin/bash
#SBATCH --partition=lotterhos
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --time=48:00:00
#SBATCH --job-name=stacks_birostris
#SBATCH --output=plink_birostris_%j.log
#SBATCH --error=plink_birostris_%j.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=eppley.m@northeastern.edu

# LD pruning for M. birostris using PLINK
# Following Humble 2023: 50 SNP window, 5 SNP step, VIF threshold of 2

module load vcftools # this is just globally available on explorer, slay

PLINK="/projects/gatins/programs_explorer/plink/plink"
VCF="minDP10_maxmiss0.8_filtInd.recode.vcf"
PREFIX="birostris_minDP10_maxmiss0.8_filtInd"

# Convert VCF to PLINK using vcftools
vcftools --vcf ${VCF} --plink --out ${PREFIX}

# Add Locus_ prefix to chromosome names in .map file
awk 'BEGIN { OFS = "\t" } { $1 = "Locus_" $1; print }' ${PREFIX}.map > ${PREFIX}.map.tmp
mv ${PREFIX}.map.tmp ${PREFIX}.map

# LD pruning
${PLINK} --file ${PREFIX} \
  --indep 50 5 2 \
  --nonfounders \
  --out ${PREFIX} \
  --allow-extra-chr \
  --debug

# Create binary PLINK files with pruned SNPs
${PLINK} --file ${PREFIX} \
  --extract ${PREFIX}.prune.in \
  --make-bed \
  --out ${PREFIX}_LDpruned \
  --allow-extra-chr

# Convert back to VCF
${PLINK} --bfile ${PREFIX}_LDpruned \
  --recode vcf \
  --out ${PREFIX}_LDpruned \
  --allow-extra-chr

echo "Total SNPs: $(wc -l < ${PREFIX}.map)"
echo "Pruned SNPs: $(wc -l < ${PREFIX}.prune.in)"
echo "Removed SNPs: $(wc -l < ${PREFIX}.prune.out)"

```
