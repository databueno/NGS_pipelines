#!/bin/bash

# identify the transcripts that correspond to the RNAseq reads in X.fastq

mkdir mapping
cd mapping
wget http://www.bioinformatics.nl/courses/BIF-30806/X.fastq
ln -s /local/data/course/genomes/Arabidopsis_thaliana/* ./
bowtie2 -x Bowtie2Index/genome -U X.fastq -S bowtie2.sam
local/prog/samtools/samtools sort -o bowtie2.bam bowtie2.sam
local/prog/samtools/samtools index bowtie2.bam

local/prog/hisat2/hisat2 -x TAIR -U X.fastq -S hisat2.sam
local/prog/samtools/samtools sort -o hisat2.bam hisat2.sam
local/prog/samtools/samtools index hisat2.bam

/local/prog/stringtie/stringtie -G genes.gtf -o X.gtf -l X hisat2.bam

/local/prog/trinity/Trinity -seqType fq --max_memory 10G -single X.fastq

makeblastdb -in genome.fa -dbtype nucl -out TAIR10
blastn -query trinity_out_dir/Trinity.fasta -db TAIR10 -evalue 1e-10 -outfmt 7 -out trinity.blast
/local/prog/kallisto/kallisto quant -i TAIR.idx -o X --single -l 200 -s 20 -t 20 X.fasta

sort -k 4 X/abundance.tsv | tail -11