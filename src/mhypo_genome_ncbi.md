## downloading the M. hypostoma genome from NCBI

first thing, get off of the login node
```
srun --partition=lotterhos --nodes=1 --ntasks=1 --cpus-per-task=4 --mem=16G --time=02:00:00 --pty bash
```

```
cd /projects/gatins/2025_Mobulid/reference_genome

wget "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/GCF_963921235.1/download?include_annotation_type=GENOME_FASTA" \
    -O ncbi_dataset.zip

unzip ncbi_dataset.zip
```

```
--2026-04-02 18:10:48--  https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/GCF_963921235.1/download?include_annotation_type=GENOME_FASTA
Connecting to 10.99.0.130:3128... connected.
Proxy request sent, awaiting response... 200 OK
Length: unspecified [application/zip]
Saving to: ‘ncbi_dataset.zip’

ncbi_dataset.zip                                  [                                                                   <=>                           ]   1.02G  37.1MB/s    in 30s     

2026-04-02 18:11:17 (35.4 MB/s) - ‘ncbi_dataset.zip’ saved [1100285896]

Archive:  ncbi_dataset.zip
  inflating: README.md               
  inflating: ncbi_dataset/data/assembly_data_report.jsonl  
  inflating: ncbi_dataset/data/GCF_963921235.1/GCF_963921235.1_sMobHyp1.1_genomic.fna  
  inflating: ncbi_dataset/data/dataset_catalog.json  
  inflating: md5sum.txt              
(base) [eppley.m@d3037 reference_genome]$ ls
md5sum.txt  ncbi_dataset  ncbi_dataset.zip  README.md
(base) [eppley.m@d3037 reference_genome]$ cd ncbi_dataset/
(base) [eppley.m@d3037 ncbi_dataset]$ ls
data
(base) [eppley.m@d3037 ncbi_dataset]$ cd data
(base) [eppley.m@d3037 data]$ ls
assembly_data_report.jsonl  dataset_catalog.json  GCF_963921235.1
```
### build indexes
then I immediately built indexes while on the srun node and moved the FASTA to the reference_genome dir

```
(base) [eppley.m@d3037 data]$
cp /projects/gatins/2025_Mobulid/reference_genome/ncbi_dataset/data/GCF_963921235.1/GCF_963921235.1_sMobHyp1.1_genomic.fna \
   /projects/gatins/2025_Mobulid/reference_genome/

cd /projects/gatins/2025_Mobulid/reference_genome
module load bwa/0.7.18
module load samtools/1.21

samtools faidx GCF_963921235.1_sMobHyp1.1_genomic.fna
bwa index GCF_963921235.1_sMobHyp1.1_genomic.fna
```


