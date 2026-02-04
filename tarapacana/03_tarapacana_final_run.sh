(base) [eppley.m@explorer-01 tarapacana]$ cat tarapacana_final_run.sh
#!/bin/bash

#SBATCH --partition=lotterhos
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --time=48:00:00
#SBATCH --job-name=tarapacana_final
#SBATCH --output=tarapacana_final_%j.log
#SBATCH --error=tarapacana_final_%j.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=eppley.m@northeastern.edu

export PATH=/projects/gatins/programs_explorer/stacks_2.68/bin:$PATH #globally available on explorer

m=3
M=2
n=2

SAMPLES_DIR="/projects/gatins/2025_Mobulid_UCSC/RAD_all_combined_bycatch/trimmed90"
BASE_DIR="/projects/gatins/2025_Mobulid/tarapacana"
OUTDIR="${BASE_DIR}/final_m${m}_M${M}_n${n}"
POPMAP="${BASE_DIR}/pop_maps/tarapacana_pop_map_3"

mkdir -p ${OUTDIR}

echo "~~~ tarapacana final run ~~~"
echo "Parameters: m=${m}, M=${M}, n=${n}"
echo "Pop map: 4 populations (EAST, OFF, CTR, WST)"
cat ${POPMAP}

ALL_SAMPLES="BYC_RMB_57 BYC_RMM_30 BYC_RMO_45 BYC_RMT_04 BYC_RMT_06 BYC_RMT_07 BYC_RMT_27 BYC_RMT_28 BYC_RMT_29 BYC_RMT_46 BYC_RMT_49 BYC_RMT_59 BYCI_RMT_69 BYCI_RMT_71"

echo "~~~ Running ustacks ~~~"
id=1
for sample in $ALL_SAMPLES; do
    ustacks -f ${SAMPLES_DIR}/${sample}.fq -o ${OUTDIR} -i $id -m $m -M $M -p 32
    ((id++))
done

echo "~~~ Running cstacks ~~~"
cstacks_cmd="cstacks -o ${OUTDIR} -p 32 -n ${n}"
for sample in $ALL_SAMPLES; do
    cstacks_cmd="$cstacks_cmd -s ${OUTDIR}/${sample}"
done
eval $cstacks_cmd

echo "~~~ Running sstacks ~~~"
sstacks -P ${OUTDIR} -M ${POPMAP} -p 32

echo "~~~ Running tsv2bam ~~~"
tsv2bam -P ${OUTDIR} -M ${POPMAP} -t 32

echo "~~~ Running gstacks ~~~"
gstacks -P ${OUTDIR} -M ${POPMAP} -t 32

echo "~~~ Running populations ~~~"
populations -P ${OUTDIR} -M ${POPMAP} -r 0.8 -p 2 --min-maf 0.05 --write-single-snp --vcf --genepop --structure --fstats --hwe -t 30

echo "~~~ Starting VCFtools ~~~"
module load vcftools

INPUT_VCF="${OUTDIR}/populations.snps.vcf"

vcftools --vcf ${INPUT_VCF} --minDP 10 --recode --recode-INFO-all --out ${OUTDIR}/minDP10

vcftools --vcf ${OUTDIR}/minDP10.recode.vcf --max-missing 0.8 --recode --recode-INFO-all --out ${OUTDIR}/minDP10_maxmiss0.8

vcftools --vcf ${OUTDIR}/minDP10_maxmiss0.8.recode.vcf --missing-indv --out ${OUTDIR}/missingness

awk '$5 > 0.4 {print $1}' ${OUTDIR}/missingness.imiss > ${OUTDIR}/remove_individuals.txt

echo "inds removed"
cat ${OUTDIR}/remove_individuals.txt

vcftools --vcf ${OUTDIR}/minDP10_maxmiss0.8.recode.vcf --remove ${OUTDIR}/remove_individuals.txt --recode --recode-INFO-all --out ${OUTDIR}/minDP10_maxmiss0.8_filtInd

echo "final VCF ${OUTDIR}/minDP10_maxmiss0.8_filtInd.recode.vcf"
grep -v "^#" ${OUTDIR}/minDP10_maxmiss0.8_filtInd.recode.vcf | wc -l
