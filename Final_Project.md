# MBIOS term project 

## Importing and prepping data files

Reference Genome
```
wget https://downloads.yeastgenome.org/sequence/S288C_reference/genome_releases/S288C_reference_genome_Current_Release.tgz

tar zxvf S288C_reference_genome_Current_Release.tgz
```

Treated and Untreated SRA files
```
gunzip -d treated_sra_data.fastq.gz

gunzip -duntreated_sra_data.fastq.gz
```
