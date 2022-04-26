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

gunzip -d untreated_data.fastq.gz
```

## Check read quality - FASTQC
```
module load fastqc 

fastqc treated_sra_data.fastq
fastqc untreated_data.fastq

#move html output files to computer via FileZilla to view in web brower
```
## Generating alignment files - bowtie2 
```
module load bowtie2

bowtie2-build -f /data/kelley/projects/eelpout/mbios/S288C_reference_genome_R64-3-1_20210421/S288C_reference_sequence_R64-3-1_20210421.fsa refgen

bowtie2 -x refgen -p 2 -U treated_sra_data.fastq -S treated_sra_data.sam

bowtie2 -x refgen -p 2 -U untreated_data.fastq.gz -S untreated_data.sam

```

## Processing alignment files- samtools
```
module load samtools 

#sam -> bam
samtools view -b treated_sra_data.sam -o treated_sra_data.bam
samtools view -b untreated_data.sam -o untreated_data.bam

#bam --> sorted.bam
samtools sort treated_sra_data.bam > treated_sra_data.sorted.bam
samtools sort untreated_data.bam > untreated_data.sorted.bam


#create pileup files 

samtools mpileup -A -f S288C_reference_genome_R64-3-1_20210421/S288C_reference_sequence_R64-3-1_20210421.fsa treated_sra_data.sorted.bam > treated_sra_data.mpileup

samtools mpileup -A -f S288C_reference_genome_R64-3-1_20210421/S288C_reference_sequence_R64-3-1_20210421.fsa untreated_data.sorted.bam > untreated_data.mpileup
```

## Scan for Variants - Varscan
```
java -jar "VarScan.v2.3.9.jar"

java -jar VarScan.v2.3.9.jar mpileup2snp treated_sra_data.mpileup --min-coverage 10 min-var-freq 0.45 --p-value 0.05 --minfreq- for-hom 0.9 > treated_varscan.txt

java -jar VarScan.v2.3.9.jar mpileup2snp untreated_data.mpileup --min-coverage 10 min-var-freq 0.45 --p-value 0.05 --min-freq- for-hom 0.9 > untreated_varscan.txt

#make files tab separated by column
paste - - < treated_varscan.txt
paste - - < untreated_varscan.txt
```
# Bedtools 
```
module load bedtools 

bedtools intersect -v -a treated_varscan_BED.txt -b untreated_varscan_BED.txt > merged.txt
bedtools intersect -wb -a merged.txt -b MBIOS_gene_positions.txt > treated_genes.txt
```
