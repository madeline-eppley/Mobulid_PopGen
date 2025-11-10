# Tarapacana RAD analysis

we are working with bycatch samples from the high seas for the species Mobula tarapacana. We have 15 RAD-sequenced samples that were collected between ~2018-2022. 

first things first! here's the endangered Mobula tarapacana, or chilean devil ray. wide-randing tropical and temperate sea species.

<img width="500" height="412" alt="image" src="https://github.com/user-attachments/assets/93c78dfe-9f71-4f1d-b10d-b7eeeed9556a" />

### ok now for the fun stuff!

#### Stacks denovo map
files uploaded from melissa and located here:
```bash
/projects/gatins/2025_Mobulid_UCSC/RAD_all_combine_bycatch
```

then we have subfolders for each species, including tarapacana

so this species is unique because we have all high-seas samples. we have lat/long info for 14/15 samples (BYC_RMT_59 is missing)

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

Let's look at some of the lat/long info from the excel sheet:


| Sample | Latitude | Longitude |
|-----------|----------|-----------|
| BYC_RMM_30 | 2°43'00.0"N | 91°59'00.0"W |
| BYC_RMB_57 | 1°26'00.0"N | 82°26'00.0"W |
| BYC-RMB-56 | 14°03'00.0"N | 82°34'00.0"W |
| BYC_RMO_45 | 2°13'00.0"N | 85°19'00.0"W |
| BYC_RMT_04 | 3°49'12.0"N | 110°04'48.0"W |
| BYC_RMT_06 | 4°16'59.9"N | 110°04'59.9"W |
| BYC_RMT_07 | 0°54'00.0"S | 92°13'48.0"W |
| BYC_RMT_27 | 2°39'00.0"N | 84°37'00.1"W |
| BYC_RMT_69 | 4°19'00.1"N | 144°24'00.0"W |
| BYC_RMT_71 | 4°36'00.0"S | 179°52'48.0"E |
| BYC_RMT_49 | 2°24'00.0"N | 90°04'48.0"W |
| BYC_RMT_28 | 3°21'00.0"N | 92°06'00.0"W |
| BYC_RMT_29 | 0°21'00.0"N | 82°19'58.8"W |
| BYC_RMT_46 | 1°46'12.0"N | 84°48'00.0"W |


<img width="1447" height="572" alt="Screenshot 2025-11-10 at 10 09 19 AM" src="https://github.com/user-attachments/assets/44fcd6ce-b5f7-419b-9e82-4acdfe730970" />


based on this output, I think a good starting place would be 3 for populations, so I made a new pop map file

this file has an eastern pacific (EST), central pacific (CTR), and western pacific group structure (WST). we are missing for BYC_RMT_59, but I assigned it eastern to start and we will go from there. 

```bash
cat > tarapacana_pop_map_3 << 'EOL'
BYC_RMB_56 EST
BYC_RMB_57 EST
BYC_RMM_30 EST
BYC_RMO_45 EST
BYC_RMT_04 CTR
BYC_RMT_06 CTR
BYC_RMT_07 EST
BYC_RMT_27 EST
BYC_RMT_28 EST
BYC_RMT_29 EST
BYC_RMT_46 EST
BYC_RMT_49 EST
BYC_RMT_59 EST
BYCI_RMT_69 CTR
BYCI_RMT_71 WST
EOL
```

I want to get an initial version of the denovo map up and working, and then we can go from there with parameter testing. 

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
