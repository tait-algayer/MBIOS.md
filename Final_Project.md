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

## Check read quality - FASTQC

module load fastqc 

fastqc treated_sra_data.fastq

# Generating alignment files - bowtie2 
```
module load bowtie2

bowtie2-build -f /data/kelley/projects/eelpout/mbios/S288C_reference_genome_R64-3-1_20210421/S288C_reference_sequence_R64-3-1_20210421.fsa refgen

bowtie2 -x refgen -p 2 -U treated_sra_data.fastq -S treated_sra_data.sam

bowtie2 -x refgen -p 2 -U untreated_sra_data.fastq.gz -S untreated_sra_data.sam
```

# Processing alignment files- samtools
```
module load samtools 

#sam -> bam
samtools view -b treated_sra_data.sam -o treated_sra_data.bam
samtools view -b untreated_sra_data.sam -o untreated_sra_data.sam

#bam --> bam.sorted
samtools sort Untreatedtrimmed.bam > Untreatedtrimmedsorted.bam
samtools sort Treatedtrimmed.bam > Treatedtrimmedsorted.bam

#create pileup files 
```
