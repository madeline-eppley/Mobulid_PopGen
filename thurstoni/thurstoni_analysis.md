# Thurstoni RAD seq analysis


### start with an estimate of the number of reads - this will give me a solid idea about the coverage before running the param opt script
I'm just going to make a quick script to get an idea about the read counts for all 36 of the thurstoni samples.

```bash
cat > /projects/gatins/2025_Mobulid/count_reads.sh << 'EOF'
#!/bin/bash
#SBATCH --partition=lotterhos
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=1:00:00
#SBATCH --job-name=count_reads
#SBATCH --output=count_reads_%j.log

# Thurstoni
echo ""
echo "=== THURSTONI ==="
echo "Sample ReadCount"
for sample in $(awk '{print $1}' /projects/gatins/2025_Mobulid_UCSC/RAD_all_combined_bycatch/pop_map_thurstoni); do
    reads=$(wc -l < /projects/gatins/2025_Mobulid_UCSC/RAD_all_combined_bycatch/trimmed90/${sample}.fq)
    reads=$((reads / 4))
    echo "${sample} ${reads}"
done | sort -k2 -n
EOF
```
and the results:
```
=== THURSTONI ===
Sample ReadCount
NIC_6003_A 481204
NIC_6016_3 636087
NIC_6003_B 751145
NIC_6003_D 783444
MT_06_ELP 787553
NIC_6016_2 787898
mthur_447 1031489
MT_05_ELP 1089071
MT_03_EC 1173626
MT_04_ELP 1242289
NIC_6003_F 1247302
MT_02_ELP 1304104
mthur_446 1306096
MT_09_EC_conc 1341550
NIC_6016_1 1401590
MT_08_EC_conc 1462728
BYC_RMO_40 1500567
NIC_2879 1547904
BYC_DESC_62 1559834
BYC_RMU_18 1605281
MT_07_ELP 1641982
MT_10_EC 1930913
MT_07_EC_conc 2073162
MT_04_EC_conc 2300461
BYC_DESC_60 2381892
NIC_6048 2525706
MT_05_EC 2606892
IN_12_MB_B 2665451
MT_01_EC 2784911
BYC_RMU_16 3201866
BYC_RMU_17 3582374
BYC_RMU_26 3819665
NIC_6003_E 4561518
BYC_RMO_48 4800917
IN_10_MB_B 5944601
BYC_RMO_47 7746021
```
it looks like 6 samples have pretty low read counts, so I'm going to remove those for building the catalog and then we can add them back in (but vcftools might filter them out later anyway)

#### samples with low read counts
NIC_6003_A, NIC_6016_3, NIC_6003_B, NIC_6003_D, MT_06_ELP, NIC_6016_2
unfortunately we do have a few from the NIC sampling location here, but we will see how they do with filtering. 

```bash
#!/bin/bash

#SBATCH --partition=lotterhos
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --time=48:00:00
#SBATCH --job-name=stacks_thurstoni
#SBATCH --output=stacks_thurstoni_%j.log
#SBATCH --error=stacks_thurstoni_%j.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=eppley.m@northeastern.edu

export PATH=/projects/gatins/programs_explorer/stacks_2.68/bin:$PATH

SAMPLES_DIR="/projects/gatins/2025_Mobulid_UCSC/RAD_all_combined_bycatch/trimmed90"
BASE_DIR="/projects/gatins/2025_Mobulid/thurstoni"
POPMAP="/projects/gatins/2025_Mobulid_UCSC/RAD_all_combined_bycatch/pop_map_thurstoni"
POPMAP_OPT="/projects/gatins/2025_Mobulid/thurstoni/pop_map_thurstoni_opt"

echo "Sample list:"
cat ${POPMAP}

# Create single-population map for optimization 
awk '{print $1 "\t" "OPT"}' ${POPMAP} > ${POPMAP_OPT}

# High coverage samples (>1M reads) for catalog building, this is excluding 6 low-coverage samples
HIGH_COV="mthur_446 mthur_447 IN_10_MB_B IN_12_MB_B MT_07_ELP MT_04_ELP MT_02_ELP MT_05_ELP NIC_6003_E NIC_6003_F NIC_6048 NIC_2879 NIC_6016_1 MT_01_EC MT_05_EC MT_04_EC_conc MT_07_EC_conc MT_10_EC MT_09_EC_conc MT_08_EC_conc MT_03_EC BYC_RMO_40 BYC_RMO_47 BYC_RMO_48 BYC_RMU_16 BYC_RMU_17 BYC_RMU_18 BYC_RMU_26 BYC_DESC_60 BYC_DESC_62"

# All samples
ALL_SAMPLES="mthur_446 mthur_447 IN_10_MB_B IN_12_MB_B MT_06_ELP MT_07_ELP MT_04_ELP MT_02_ELP MT_05_ELP NIC_6003_B NIC_6003_D NIC_6003_A NIC_6016_2 NIC_6016_3 NIC_6003_E NIC_6048 NIC_2879 NIC_6003_F NIC_6016_1 MT_01_EC MT_05_EC MT_04_EC_conc MT_07_EC_conc MT_10_EC MT_09_EC_conc MT_08_EC_conc MT_03_EC BYC_RMO_40 BYC_RMO_47 BYC_RMO_48 BYC_RMU_16 BYC_RMU_17 BYC_RMU_18 BYC_RMU_26 BYC_DESC_60 BYC_DESC_62"

mkdir -p ${BASE_DIR}/opt

m=3

# Paris et al. 2017: M=2-5 is testing range
# Testing M=2,3,4,5 with n=M-1, n=M, n=M+1
for M in 2 3 4 5; do
    
    echo "=== Running ustacks with m=${m}, M=${M} ==="
    USTACKS_DIR="${BASE_DIR}/opt/ustacks_m${m}_M${M}"
    mkdir -p ${USTACKS_DIR}
    
    id=1
    for sample in $ALL_SAMPLES; do
        ustacks -f ${SAMPLES_DIR}/${sample}.fq -o ${USTACKS_DIR} -i $id -m $m -M $M -p 32
        ((id++))
    done
    
    for n in $((M - 1)) $M $((M + 1)); do
        
        echo "=== Building catalog and running populations: m=${m}, M=${M}, n=${n} ==="
        OUT_DIR="${BASE_DIR}/opt/m${m}_M${M}_n${n}"
        mkdir -p ${OUT_DIR}
        
        # copy ustacks output
        cp ${USTACKS_DIR}/*.tsv ${OUT_DIR}/
        
        # cstacks - HIGH COVERAGE samples only
        cstacks_cmd="cstacks -o ${OUT_DIR} -p 32 -n ${n}"
        for sample in $HIGH_COV; do
            cstacks_cmd="$cstacks_cmd -s ${OUT_DIR}/${sample}"
        done
        eval $cstacks_cmd
        
        # sstacks - ALL samples
        for sample in $ALL_SAMPLES; do
            sstacks -c ${OUT_DIR}/catalog -s ${OUT_DIR}/${sample} -o ${OUT_DIR} -p 32
        done
        
        # tsv2bam and gstacks
        tsv2bam -P ${OUT_DIR} -M ${POPMAP_OPT} -t 32
        gstacks -P ${OUT_DIR} -M ${POPMAP_OPT} -t 32
        
        # populations with r80 for optimization (single pop, -p 1)
        populations -P ${OUT_DIR} -M ${POPMAP_OPT} -r 0.8 -p 1 -t 30
        
    done
done

ls -l ${BASE_DIR}/opt/
```
