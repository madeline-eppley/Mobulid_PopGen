# Tarapacana RADseq analysis

we are working with bycatch samples from the high seas for the species Mobula tarapacana. We have 15 RAD-sequenced samples that were collected between ~2018-2022. 

first things first! here's the [critically endangered](https://www.iucnredlist.org/es/species/60199/279077133) Mobula tarapacana, or the sicklefin devil ray. wide-ranging tropical/temperate species.

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

## optimizing parameters
ok time to get serious. we've got stacks up and running and now we should focus on the parameter optimization. We can determine the optimal params from the highest yield loci at the end of the pipeline. 

#### stacks params
**Options:**
`-m` — Minimum depth of coverage required to create a stack (default 3).  
`-M` — number of mismatches allowed between stacks within individuals (for ustacks).  
`-n` — number of mismatches allowed between stacks between individuals (for cstacks).  
`-samples [path]` — specify a path to the directory of samples (samples will be read from population map).  
`--popmap [path]` — path to a population map file (format is "[name] TAB [pop]", one sample per line).  
`-o [path]` — path to write pipeline output files.  
`-p`,`--min-populations` — minimum number of populations a locus must be present in to process a locus (for populations; default: 1). 
`-r`,`--min-samples-per-pop` — minimum percentage of individuals in a population required to process a locus for that population (for populations; default: 0). 

#### Remy and I agree that -m 3 and -M 3 are default, let's keep those and instead vary parameter -n
for some further proof beyond the BSB dataset, here's what [Paris et al 2017](https://besjournals.onlinelibrary.wiley.com/doi/epdf/10.1111/2041-210X.12775) says about stacks parameter space across 4 different datasets: "We recommend setting m high enough only to deny errantreads the status of a putative allele. We therefore do not sug-gest using a value for m > 5 and we have demonstrated herethat the default value of m3 was favourable for all test datasets"

tldr: -m 3 is default and should be good here

```bash
(base) [eppley.m@explorer-02 tarapacana]$ cat tarapacana_opt_n.sh 
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

# new dir for our opt folder
mkdir -p /projects/gatins/2025_Mobulid/tarapacana/opt

# parameters to test for n
params="
2
3
4
5
6
7"

# now denovo_map.pl for each n value
for p in $params
do
    echo "Running with n = ${p}"
    mkdir -p /projects/gatins/2025_Mobulid/tarapacana/opt/n${p}
    
    denovo_map.pl \
        -m 3 \
        -M 3 \
        -n ${p} \
        -T 32 \
        -o /projects/gatins/2025_Mobulid/tarapacana/opt/n${p} \
        --popmap /projects/gatins/2025_Mobulid_UCSC/RAD_all_combined_bycatch/pop_map_tarapacana \
        --samples /projects/gatins/2025_Mobulid_UCSC/RAD_all_combined_bycatch/trimmed90
done
```



## testing values of p (p = pops)
now we want to get into population testing since we don't have pre-existing expectations about # of pops in the data set. 
for our optimal catalog with `-m 3 -M 3 -n opt`, we are going to cycle through the populations with different pop maps that each have different groups of samples based on collection location. 

#### let's organize 🫡
step 1: move a copy of the p=1 og pop map into the new pop map folder and rename it to be pop_map_1 to make my life easier

`cp /projects/gatins/2025_Mobulid_UCSC/RAD_all_combined_bycatch/pop_map_tarapacana /projects/gatins/2025_Mobulid/tarapacana/pop_maps/tarapacana_pop_map_1`

and now repeat to do the same with our new n=3 pop map that we made earlier based on the high seas sampling locations

`mv /projects/gatins/2025_Mobulid/tarapacana/tarapacana_pop_map_3 /projects/gatins/2025_Mobulid/tarapacana/pop_maps/`

let's make a pop map file for p=2 by combining central with western:

```bash
cat > /projects/gatins/2025_Mobulid/tarapacana/pop_maps/tarapacana_pop_map_2 << 'EOF'
BYC_RMB_56 EST
BYC_RMB_57 EST
BYC_RMM_30 EST
BYC_RMO_45 EST
BYC_RMT_04 WST
BYC_RMT_06 WST
BYC_RMT_07 EST
BYC_RMT_27 EST
BYC_RMT_28 EST
BYC_RMT_29 EST
BYC_RMT_46 EST
BYC_RMT_49 EST
BYC_RMT_59 EST
BYCI_RMT_69 WST
BYCI_RMT_71 WST
EOF
```

so at this point, p=2 and p=3 are kind of the end of my educated "guessing" based on geography. we will also examine the PCA with p=1 to check for signs of population structure. the last thing we can do is make every individual their OWN pop, and check the PCA there too. so let's make that pop map where every ind is its own pop:

```bash
cat > /projects/gatins/2025_Mobulid/tarapacana/pop_maps/pop_map_15 << 'EOF'
BYCI_RMT_69 1
BYCI_RMT_71 2
BYC_RMB_56 3
BYC_RMB_57 4
BYC_RMM_30 5
BYC_RMO_45 6
BYC_RMT_04 7
BYC_RMT_06 8
BYC_RMT_07 9
BYC_RMT_27 10
BYC_RMT_28 11
BYC_RMT_29 12
BYC_RMT_46 13
BYC_RMT_49 14
BYC_RMT_59 15
EOF
```


### running populations for each pop map
```bash
#!/bin/bash

#SBATCH --partition=lotterhos
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --time=48:00:00
#SBATCH --job-name=tarapacana_pops
#SBATCH --output=tarapacana_pops_%j.log
#SBATCH --error=tarapacana_pops_%j.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=eppley.m@northeastern.edu

export PATH=/projects/gatins/programs_explorer/stacks_2.68/bin:$PATH

OPT_N=X

mkdir -p /projects/gatins/2025_Mobulid/tarapacana/populations

populations \
  -P /projects/gatins/2025_Mobulid/tarapacana/opt/n${OPT_N} \
  -M /projects/gatins/2025_Mobulid/tarapacana/pop_maps/pop_map_15 \
  -r 0.80 \
  -p 15 \
  --min-maf 0.05 \
  --write-single-snp \
  --vcf \
  --genepop \
  --structure \
  --fstats \
  --hwe \
  -O /projects/gatins/2025_Mobulid/tarapacana/populations/15_pops \
  -t 30

populations \
  -P /projects/gatins/2025_Mobulid/tarapacana/opt/n${OPT_N} \
  -M /projects/gatins/2025_Mobulid/tarapacana/pop_maps/tarapacana_pop_map_1 \
  -r 0.80 \
  -p 1 \
  --min-maf 0.05 \
  --write-single-snp \
  --vcf \
  --genepop \
  --structure \
  --fstats \
  --hwe \
  -O /projects/gatins/2025_Mobulid/tarapacana/populations/1_pops \
  -t 30

populations \
  -P /projects/gatins/2025_Mobulid/tarapacana/opt/n${OPT_N} \
  -M /projects/gatins/2025_Mobulid/tarapacana/pop_maps/tarapacana_pop_map_2 \
  -r 0.80 \
  -p 2 \
  --min-maf 0.05 \
  --write-single-snp \
  --vcf \
  --genepop \
  --structure \
  --fstats \
  --hwe \
  -O /projects/gatins/2025_Mobulid/tarapacana/populations/2_pops \
  -t 30

populations \
  -P /projects/gatins/2025_Mobulid/tarapacana/opt/n${OPT_N} \
  -M /projects/gatins/2025_Mobulid/tarapacana/pop_maps/tarapacana_pop_map_3 \
  -r 0.80 \
  -p 3 \
  --min-maf 0.05 \
  --write-single-snp \
  --vcf \
  --genepop \
  --structure \
  --fstats \
  --hwe \
  -O /projects/gatins/2025_Mobulid/tarapacana/populations/3_pops \
  -t 30
```

## pop structure time!!!


