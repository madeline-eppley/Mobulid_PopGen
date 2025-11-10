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



SAMPLE_COUNT=$(wc -l < /projects/gatins/2025_Mobulid_UCSC/RAD_all_combined_bycatch/pop_map_tarapacana)
echo "Total num of samples: $SAMPLE_COUNT"

echo "Sample list:"
cat /projects/gatins/2025_Mobulid_UCSC/RAD_all_combined_bycatch/pop_map_tarapacana

echo "Starting denovo_map.pl"

denovo_map.pl \
  -m 3 \
  -M 3 \
  -n 4 \
  -T 32 \
  -o /projects/gatins/2025_Mobulid/tarapacana \
  --popmap /projects/gatins/2025_Mobulid_UCSC/RAD_all_combined_bycatch/pop_map_tarapacana \
  --samples /projects/gatins/2025_Mobulid_UCSC/RAD_all_combined_bycatch/trimmed90



ls -l /projects/gatins/2025_Mobulid/tarapacana
