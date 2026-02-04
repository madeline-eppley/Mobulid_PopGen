(base) [eppley.m@explorer-01 tarapacana]$ cat run_populations_gal.sh
#!/bin/bash
#SBATCH --partition=short
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --job-name=tar_pop_gal
#SBATCH --output=tar_pop_gal_%j.log
#SBATCH --error=tar_pop_gal_%j.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=eppley.m@northeastern.edu

export PATH=/projects/gatins/programs_explorer/stacks_2.68/bin:$PATH # globally available on explorer

BASE_DIR="/projects/gatins/2025_Mobulid/tarapacana"
STACKS_DIR="${BASE_DIR}/final_m3_M2_n2"
OUTDIR="${BASE_DIR}/populations_gal"
POPMAP="${BASE_DIR}/popmap_gal.txt"

mkdir -p ${OUTDIR}

echo "~~~ Running populations with new 4 pop (EAST, GAL, CTR, WST) popmap ~~~"
echo "Input: ${STACKS_DIR}"
echo "Output: ${OUTDIR}"
echo "Popmap:"
cat ${POPMAP}

populations -P ${STACKS_DIR} \
  -M ${POPMAP} \
  -O ${OUTDIR} \
  -r 0.8 \
  -p 2 \
  --min-maf 0.05 \
  --write-single-snp \
  --vcf \
  --fstats \
  -t 8

echo "~~~ Populations complete, starting VCFtools filtering ~~~"

module load vcftools

INPUT_VCF="${OUTDIR}/populations.snps.vcf"

# minDP 10
vcftools --vcf ${INPUT_VCF} --minDP 10 --recode --recode-INFO-all --out ${OUTDIR}/minDP10

# max-missing 0.8
vcftools --vcf ${OUTDIR}/minDP10.recode.vcf --max-missing 0.8 --recode --recode-INFO-all --out ${OUTDIR}/minDP10_maxmiss0.8

# check missingness per ind
vcftools --vcf ${OUTDIR}/minDP10_maxmiss0.8.recode.vcf --missing-indv --out ${OUTDIR}/missingness

echo "~~~ Individual missingness ~~~"
cat ${OUTDIR}/missingness.imiss

# remove individuals with >40% missing
awk '$5 > 0.4 {print $1}' ${OUTDIR}/missingness.imiss > ${OUTDIR}/remove_individuals.txt

if [ -s ${OUTDIR}/remove_individuals.txt ]; then
    echo "Removing individuals with >40% missing data:"
    cat ${OUTDIR}/remove_individuals.txt
    vcftools --vcf ${OUTDIR}/minDP10_maxmiss0.8.recode.vcf --remove ${OUTDIR}/remove_individuals.txt --recode --recode-INFO-all --out ${OUTDIR}/tarapacana_4pop_final
else
    echo "No individuals to remove, copying final VCF"
    cp ${OUTDIR}/minDP10_maxmiss0.8.recode.vcf ${OUTDIR}/tarapacana_4pop_final.recode.vcf
fi

echo "final SNP count"
grep -v "^#" ${OUTDIR}/tarapacana_4pop_final.recode.vcf | wc -l

ls -la ${OUTDIR}
