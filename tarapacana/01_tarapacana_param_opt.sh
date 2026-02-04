```bash
(base) [eppley.m@explorer-01 tarapacana]$ cat tarapacana_param_opt.sh 
#!/bin/bash

#SBATCH --partition=lotterhos
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --time=48:00:00
#SBATCH --job-name=stacks_tarapacana
#SBATCH --output=stacks_tarapacana_%j.log
#SBATCH --error=stacks_tarapacana_%j.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=eppley.m@northeastern.edu

export PATH=/projects/gatins/programs_explorer/stacks_2.68/bin:$PATH

SAMPLES_DIR="/projects/gatins/2025_Mobulid_UCSC/RAD_all_combined_bycatch/trimmed90"
BASE_DIR="/projects/gatins/2025_Mobulid/tarapacana"
POPMAP="/projects/gatins/2025_Mobulid/tarapacana/pop_map_tarapacana_nodup"
POPMAP_OPT="/projects/gatins/2025_Mobulid/tarapacana/pop_map_tarapacana_opt"

# BYC_RMB_56 removed - duplicate/highly related to BYC_RMB_57 (kinship ~0.498)
# BYC_RMB_56 had higher missingness (~4% vs 0.7%)
# all remaining samples are >10x coverage, so no high-coverage filtering needed
cat ${POPMAP}

# single pop map
awk '{print $1 "\t" "OPT"}' ${POPMAP} > ${POPMAP_OPT}

# all samples
ALL_SAMPLES="BYC_RMB_57 BYC_RMM_30 BYC_RMO_45 BYC_RMT_04 BYC_RMT_06 BYC_RMT_07 BYC_RMT_27 BYC_RMT_28 BYC_RMT_29 BYC_RMT_46 BYC_RMT_49 BYC_RMT_59 BYCI_RMT_69 BYCI_RMT_71"

mkdir -p ${BASE_DIR}/opt

m=3

# Paris et al. 2017: M=2-5 is testing range
# testing M=2,3,4,5 with n=M-1, n=M, n=M+1
for M in 2 3 4 5; do
    
    echo "~~~ Running ustacks with m=${m}, M=${M} ~~~"
    USTACKS_DIR="${BASE_DIR}/opt/ustacks_m${m}_M${M}"
    mkdir -p ${USTACKS_DIR}
    
    id=1
    for sample in $ALL_SAMPLES; do
        ustacks -f ${SAMPLES_DIR}/${sample}.fq -o ${USTACKS_DIR} -i $id -m $m -M $M -p 32
        ((id++))
    done
    
    for n in $((M - 1)) $M $((M + 1)); do
        
        echo "~~~ Building catalog and running populations: m=${m}, M=${M}, n=${n} ~~~"
        OUT_DIR="${BASE_DIR}/opt/m${m}_M${M}_n${n}"
        mkdir -p ${OUT_DIR}
        
        # copy ustacks output
        cp ${USTACKS_DIR}/*.tsv ${OUT_DIR}/
        
        # cstacks - ALL samples (all >10x coverage)
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
