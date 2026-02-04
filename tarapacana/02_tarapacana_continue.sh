```bash
#!/bin/bash

#SBATCH --partition=lotterhos
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --time=24:00:00
#SBATCH --job-name=tarapacana_cont
#SBATCH --output=tarapacana_cont_%j.log
#SBATCH --error=tarapacana_cont_%j.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=eppley.m@northeastern.edu

export PATH=/projects/gatins/programs_explorer/stacks_2.68/bin:$PATH

BASE_DIR="/projects/gatins/2025_Mobulid/tarapacana"
POPMAP_OPT="/projects/gatins/2025_Mobulid/tarapacana/pop_map_tarapacana_opt"

# loop through existing parameter directories to check if completed in case of time out
for OUT_DIR in ${BASE_DIR}/opt/m3_M*_n*; do
    echo "~~~ start ${OUT_DIR} ~~~"
    
    # sstacks - using correct -P syntax
    sstacks -P ${OUT_DIR} -M ${POPMAP_OPT} -p 32
    
    # tsv2bam
    tsv2bam -P ${OUT_DIR} -M ${POPMAP_OPT} -t 32
    
    # gstacks
    gstacks -P ${OUT_DIR} -M ${POPMAP_OPT} -t 32
    
    # populations
    populations -P ${OUT_DIR} -M ${POPMAP_OPT} -r 0.8 -p 1 -t 30
    
    echo "~~~ done with ${OUT_DIR} ~~~"
done
```
