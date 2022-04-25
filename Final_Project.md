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
```
module load fastqc 

fastqc treated_sra_data.fastq
fastqc untreated_sra_data.fastq

#move html output files to computer via FileZilla to view in web brower
```
## Generating alignment files - bowtie2 
```
module load bowtie2

bowtie2-build -f /data/kelley/projects/eelpout/mbios/S288C_reference_genome_R64-3-1_20210421/S288C_reference_sequence_R64-3-1_20210421.fsa refgen

bowtie2 -x refgen -p 2 -U treated_sra_data.fastq -S treated_sra_data.sam

bowtie2 -x refgen -p 2 -U untreated_sra_data.fastq.gz -S untreated_sra_data.sam

bowtie2 -x refgen -p 2 -U sra_data.fastq -S sra_data.fastq2.sam
```

## Processing alignment files- samtools
```
module load samtools 

#sam -> bam
samtools view -b treated_sra_data.sam -o treated_sra_data.bam
samtools view -b untreated_sra_data.sam -o untreated_sra_data.sam

samtools view -b 2sra_data.sam -o 2sra_data.sam

#bam --> bam.sorted
samtools sort Untreatedtrimmed.bam > Untreatedtrimmedsorted.bam
samtools sort Treatedtrimmed.bam > Treatedtrimmedsorted.bam

#create pileup files 

samtools mpileup -A -f S288C_reference_genome_R64-3-1_20210421/S288C_reference_sequence_R64-3-1_20210421.fsa treated_sra_data.bam > treated_sra_data.mpileup

samtools mpileup -A -f S288C_reference_genome_R64-3-1_20210421/S288C_reference_sequence_R64-3-1_20210421.fsa untreated_sra_data.bam > untreated_sra_data.mpileup
```

## Scan for Variants - Varscan
```
java -jar "VarScan.v2.3.9.jar"

java -jar VarScan.jar mpileup2snp Untreated.mpileup --min-coverage 10 min-var-freq 0.45 --p-value 0.05 --minfreq-
for-hom 0.9 > UntreatedS288C_SNPs_Varscan.txt
java -jar VarScan.jar mpileup2snp Treated2.mpileup --min-coverage 10 min-var-freq 0.45 --p-value 0.05 --min-freqfor-
hom 0.9 > Treated2S288C_SNPs_Varscan.txt
```
# Bedtools 
```
bedtools intersect -v -a Treated2S288C_SNPs_Varscan.txt -b UntreatedS288C_SNPs_Varscan.txt >
Uniquetreated.txt
bedtools intersect -wb -a Uniquetreated.txt -b augustusGene.txt > Treateduniquegenes.txt
```
