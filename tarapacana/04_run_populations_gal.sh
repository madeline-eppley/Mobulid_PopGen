(base) [eppley.m@explorer-01 tarapacana]$ cat run_populations_gal.sh
#!/bin/bash

#SBATCH --partition=lotterhos
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

export PATH=/projects/gatins/programs_explorer/stacks_2.68/bin:$PATH #globally avail on explorer

BASE_DIR="/projects/gatins/2025_Mobulid/tarapacana"
STACKS_DIR="${BASE_DIR}/final_m3_M2_n2"
OUTDIR="${BASE_DIR}/populations_gal"
POPMAP="${BASE_DIR}/popmap_gal.txt"

mkdir -p ${OUTDIR}

echo "~~~ Running populations with new 4 pop (EAST, Central, West, Galapagos/GAL) popmap ~~~"
echo "Input: ${STACKS_DIR}"
echo "Output: ${OUTDIR}"
echo "Popmap:"
cat ${POPMAP}

populations -P ${STACKS_DIR} \
  -M ${POPMAP} \
  -O ${OUTDIR} \
  -r 0.8 \
  --min-maf 0.05 \
  --vcf \
  --fstats \
  -t 8

ls -la ${OUTDIR}
