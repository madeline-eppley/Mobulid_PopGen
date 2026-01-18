# Mobular RAD seq analysis

### quick ballpark on read counts

```
=== MOBULAR ===
Sample ReadCount
BYC_RMM_23 1569296
BYC_RMM_58 1660809
BYC_RMU_42 1935879
BYC_DESC_55 2091939
BYC_DESC_63 2452890
BYC_RMM_31 2502172
BYC_RMM_22 2783324
BYCI_RMM_66 3148974
BYC_RMM_51 3251828
BYC_RMM_44 3366970
BYC_RMM_14 3459831
BYC_RMU_41 3517749
BYC_RMM_05 3982051
BYCI_DESC_67 4021644
BYCI_RMM_68 4054104
BYC_RMM_35 4081332
BYC_DESC_54 4707895
BYC_RMM_11 4791938
BYC_RMM_13 4996620
BYC_RMM_33 5424949
BYC_RMB_34 6020354
BYC_RMM_50 6227607
BYC_RMM_38 6444491
BYCI_RMM_72 6963395
BYC_RMM_12 6998277
BYC_RMM_10 7263832
BYC_RMM_24 7667013
BYC_RMM_32 7694619
BYC_RMM_03 8069053
BYC_DESC_64 8566863
BYC_RMM_21 9396439
BYC_RMB_36 9708388
BYC_RMM_20 10451314
BYC_RMM_02 13254849
BYC_RMM_25 15347195
```


### param opt script
```bash
#!/bin/bash

#SBATCH --partition=lotterhos
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --time=48:00:00
#SBATCH --job-name=stacks_mobular
#SBATCH --output=stacks_mobular_%j.log
#SBATCH --error=stacks_mobular_%j.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=eppley.m@northeastern.edu

export PATH=/projects/gatins/programs_explorer/stacks_2.68/bin:$PATH

SAMPLES_DIR="/projects/gatins/2025_Mobulid_UCSC/RAD_all_combined_bycatch/trimmed90"
BASE_DIR="/projects/gatins/2025_Mobulid/mobular"
POPMAP="/projects/gatins/2025_Mobulid_UCSC/RAD_all_combined_bycatch/pop_map_mobular"
POPMAP_OPT="/projects/gatins/2025_Mobulid/mobular/pop_map_mobular_opt"

echo "Sample list:"
cat ${POPMAP}

# Create single-population map for optimization
awk '{print $1 "\t" "OPT"}' ${POPMAP} > ${POPMAP_OPT}

# Build catalog on all samples
ALL_SAMPLES="BYC_DESC_54 BYC_DESC_55 BYC_DESC_63 BYC_DESC_64 BYCI_DESC_67 BYC_RMB_34 BYC_RMB_36 BYC_RMM_02 BYC_RMM_03 BYC_RMM_05 BYC_RMM_10 BYC_RMM_11 BYC_RMM_12 BYC_RMM_13 BYC_RMM_14 BYC_RMM_20 BYC_RMM_21 BYC_RMM_22 BYC_RMM_23 BYC_RMM_24 BYC_RMM_25 BYC_RMM_31 BYC_RMM_32 BYC_RMM_33 BYC_RMM_35 BYC_RMM_38 BYC_RMM_44 BYC_RMM_50 BYC_RMM_51 BYC_RMM_58 BYCI_RMM_66 BYCI_RMM_68 BYCI_RMM_72 BYC_RMU_41 BYC_RMU_42"

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
        
        # cstacks - ALL samples
        cstacks_cmd="cstacks -o ${OUT_DIR} -p 32 -n ${n}"
        for sample in $ALL_SAMPLES; do
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
