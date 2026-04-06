(base) [eppley.m@explorer-02 all_mobulids]$ cat all_mobulids_run.sh 
#!/bin/bash
#SBATCH --partition=lotterhos
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=256G
#SBATCH --time=96:00:00
#SBATCH --job-name=all_mobulids
#SBATCH --output=all_mobulids_%j.log
#SBATCH --error=all_mobulids_%j.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=eppley.m@northeastern.edu

export PATH=/projects/gatins/programs_explorer/stacks_2.68/bin:$PATH

# goal is to create a comparative dataset that has all species in it with one filtering process for effective pop size, phylogenetics

# chosing pre-defined optimal parameters that were optimal across most species 
m=3
M=2
n=2

SAMPLES_DIR="/projects/gatins/2025_Mobulid_UCSC/RAD_all_combined_bycatch/trimmed90"
BASE_DIR="/projects/gatins/2025_Mobulid/all_mobulids"
OUTDIR="${BASE_DIR}/final_m${m}_M${M}_n${n}"
POPMAP="${BASE_DIR}/all_mobulids_popmap.txt"

mkdir -p ${OUTDIR}

echo "~~~ all mobulids comparative stacks run ~~~"
echo "Start time: $(date)"

# count samples per species
echo "~~~~~ pop map summary ~~~"
echo "total samples: $(wc -l < ${POPMAP})"
for sp in birostris munkiana thurstoni mobular tarapacana; do
    echo "  ${sp}: $(grep -c "${sp}" ${POPMAP})"
done
echo ""
cat ${POPMAP}
echo ""

# these are all 109 individuals that we have in the species specific datasets
# birostris (22)
BIROSTRIS="BYC_RMB_01 IN_1_MB IN_2_MB IN_4_MB IN_5_MB IN_6_MB IN_7_MB PER_001_MB PER_003_MB PER_004_MB PER_005_MB PER_006_MB PER_007_MB PER_DZW81_4_MB REV_10_MB REV_13_MB REV_14_MB_B REV_15_MB REV_17_MB REV_18_MB REV_19_MB REV_20_MB_B"

# munkiana (13)
MUNKIANA="mu_11 mu_12 mu_21 mu_3 mu_4 mu_7 NIC_2728_01 NIC_2740_01 PER_DCLCW130_MM PER_DCLCW87_MM PER_DPPW019_MM PER_DPPW28_MM PER_DPPW98_MM"

# thurstoni (28)
THURSTONI="mthur_446 IN_10_MB_B IN_12_MB_B MT_07_ELP MT_04_ELP MT_02_ELP MT_05_ELP NIC_6003_E NIC_6048 NIC_6003_F NIC_6016_1 MT_01_EC MT_05_EC MT_04_EC_conc MT_07_EC_conc MT_10_EC MT_09_EC_conc MT_08_EC_conc MT_03_EC BYC_RMO_40 BYC_RMO_47 BYC_RMO_48 BYC_RMU_16 BYC_RMU_17 BYC_RMU_18 BYC_RMU_26 BYC_DESC_60 BYC_DESC_62"

# mobular (33)
MOBULAR="BYC_RMM_10 BYC_RMM_11 BYC_RMM_12 BYC_RMM_13 BYC_RMM_14 BYC_RMM_20 BYC_RMM_21 BYC_RMM_22 BYC_RMM_24 BYC_RMM_25 BYC_RMM_31 BYC_RMM_32 BYC_RMB_34 BYC_RMM_35 BYC_RMB_36 BYC_RMM_38 BYC_RMU_41 BYC_RMU_42 BYC_RMM_44 BYC_RMM_50 BYC_RMM_51 BYC_DESC_54 BYC_DESC_55 BYC_RMM_58 BYC_DESC_63 BYC_DESC_64 BYC_RMM_02 BYC_RMM_03 BYC_RMM_05 BYCI_RMM_68 BYCI_RMM_66 BYCI_DESC_67 BYCI_RMM_72"

# tarapacana (13)
TARAPACANA="BYC_RMB_57 BYC_RMO_45 BYC_RMT_27 BYC_RMT_29 BYC_RMT_46 BYC_RMM_30 BYC_RMT_07 BYC_RMT_28 BYC_RMT_49 BYC_RMT_04 BYC_RMT_06 BYCI_RMT_69 BYCI_RMT_71"

ALL_SAMPLES="${BIROSTRIS} ${MUNKIANA} ${THURSTONI} ${MOBULAR} ${TARAPACANA}"

# ~~~~ ustacks~~~~~~
echo "~~~~ ustacks~~~~~~"
echo "ustacks start: $(date)"
id=1
for sample in $ALL_SAMPLES; do
    echo "  ustacks: ${sample} (id=${id}) -- $(date)"
    ustacks -f ${SAMPLES_DIR}/${sample}.fq -o ${OUTDIR} -i $id -m $m -M $M -p 32
    ((id++))
done
echo "~~~~ ustacks~~~~~~ complete: $(date)"
echo ""

# ~~~~ cstacks~~~~~~
echo "~~~~ cstacks~~~~~~"
echo "cstacks start: $(date)"
cstacks_cmd="cstacks -o ${OUTDIR} -p 32 -n ${n}"
for sample in $ALL_SAMPLES; do
    cstacks_cmd="$cstacks_cmd -s ${OUTDIR}/${sample}"
done
eval $cstacks_cmd
echo "cstacks complete: $(date)"
echo ""

# ~~~~ sstacks~~~~~~
echo "~~~~ sstacks~~~~~~"
echo "sstacks start: $(date)"
sstacks -P ${OUTDIR} -M ${POPMAP} -p 32
echo "sstacks complete: $(date)"
echo ""

# ~~~~ tsv2bam~~~~~~
echo "~~~~~ tsv2bam~~~~~~"
echo "tsv2bam start: $(date)"
tsv2bam -P ${OUTDIR} -M ${POPMAP} -t 32
echo "tsv2bam complete: $(date)"
echo ""

# ~~~~ gstacks ~~~~~~
echo " ~~~~ gstacks ~~~~~~"
echo "gstacks start: $(date)"
gstacks -P ${OUTDIR} -M ${POPMAP} -t 32
echo "gstacks complete: $(date)"
echo ""

#  ~~~~ pops ~~~~~~
echo " ~~~~ pops ~~~~~~"
echo "pops start: $(date)"
populations -P ${OUTDIR} \
  -M ${POPMAP} \
  -r 0.8 \
  -p 2 \
  --min-maf 0.05 \
  --write-single-snp \
  --vcf \
  --genepop \
  --structure \
  --fstats \
  --hwe \
  -t 30
echo "populations complete: $(date)"
echo ""

#  ~~~~ vcftools ~~~~~~
echo "~~~~~~~ vcftools ~~~~~~~~~"
module load vcftools

INPUT_VCF="${OUTDIR}/populations.snps.vcf"

# minDP10
vcftools --vcf ${INPUT_VCF} --minDP 10 --recode --recode-INFO-all --out ${OUTDIR}/minDP10

# max-missing 0.8
vcftools --vcf ${OUTDIR}/minDP10.recode.vcf --max-missing 0.8 --recode --recode-INFO-all --out ${OUTDIR}/minDP10_maxmiss0.8

# check missingness per ind
vcftools --vcf ${OUTDIR}/minDP10_maxmiss0.8.recode.vcf --missing-indv --out ${OUTDIR}/missingness

echo "~~~ inds missingness ~~~"
cat ${OUTDIR}/missingness.imiss

# remove inds with >40% missing
awk '$5 > 0.4 {print $1}' ${OUTDIR}/missingness.imiss > ${OUTDIR}/remove_individuals.txt

if [ -s ${OUTDIR}/remove_individuals.txt ]; then
    echo "removed individuals with >40% missing data:"
    cat ${OUTDIR}/remove_individuals.txt
    vcftools --vcf ${OUTDIR}/minDP10_maxmiss0.8.recode.vcf \
      --remove ${OUTDIR}/remove_individuals.txt \
      --recode --recode-INFO-all \
      --out ${OUTDIR}/all_mobulids_final
else
    echo "no inds with high missingness"
    cp ${OUTDIR}/minDP10_maxmiss0.8.recode.vcf ${OUTDIR}/all_mobulids_final.recode.vcf
fi

echo "snp count:"
grep -v "^#" ${OUTDIR}/all_mobulids_final.recode.vcf | wc -l
echo "inds count"
grep "^#CHROM" ${OUTDIR}/all_mobulids_final.recode.vcf | tr '\t' '\n' | tail -n +10 | wc -l
echo "end time: $(date)"
