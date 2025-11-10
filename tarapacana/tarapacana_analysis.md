# Analysis steps for high seas samples tarapacana

## Stacks denovo map
files uploaded from melissa at 
```bash
/projects/gatins/2025_Mobulid_UCSC/RAD_all_combine_bycatch
```

then we have subfolders for each species, including tarapacana

so this species is unique because we have all high-seas samples. we do have some lat/long info for some of the samples, however

`pop_map_tarapacana`

```bash
BYC_RMB_56	BYC
BYC_RMB_57	BYC
BYC_RMM_30	BYC
BYC_RMO_45	BYC
BYC_RMT_04	BYC
BYC_RMT_06	BYC
BYC_RMT_07	BYC
BYC_RMT_27	BYC
BYC_RMT_28	BYC
BYC_RMT_29	BYC
BYC_RMT_46	BYC
BYC_RMT_49	BYC
BYC_RMT_59	BYC
BYCI_RMT_69	BYC
BYCI_RMT_71	BYC
```

Let's look at some of the lat/long info from the excel sheet (full version uploaded in /tarapacana
BYC_RMM_30
2°43'00.0"N 91°59'00.0"W

BYC_RMB_57
1°26'00.0"N 82°26'00.0"W

BYC-RMB-56
14°03'00.0"N 82°34'00.0"W



Inital test run of this denovo map script works! 

```bash
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

```
