#!/bin/bash

# Modify path for your own computer/system/machine
input_dir="/media/dario/Data/HG01-prova/INPUT/fastq"
ref_dir="/media/dario/Data/HG01-prova/INPUT/hg38"
output_dir="/media/dario/Data/HG01-prova/OUTPUT"
bwa_output="/media/dario/Data/HG01-prova/INPUT/bwa_output"
fastqc_dir="/media/dario/Data/HG01-prova/INPUT/fastqc_ouput"
supporting_files_dir="/media/dario/Data/HG01-prova/INPUT/supporting_files/hg38"
data_dir="/media/dario/Data/HG01-prova/INPUT/data"

for sample in $(ls $input_dir | cut -f1 -d'.' | sort -u); do
    #producing fastqc report html  for R1 and R2 
    fastqc $input_dir/$sample.*.R1.fastq -o $fastqc_dir
    fastqc $input_dir/$sample.*.R2.fastq -o $fastqc_dir
    #alignment to the references
    bwa mem -t 8 -R "@RG\tID:${sample}\tSM:${sample}\tLB:${sample}\tPU:${sample}\tPL:${sample}" ${ref_dir}/hg38.fasta $input_dir/$sample.*.R1.fastq $input_dir/$sample.*.R2.fastq -O "$bwa_output/$sample.paired.sam"
    #mark duplicates in the alignment phase 
    gatk MarkDuplicatesSpark -I "$bwa_output/$sample.paired.sam" -O "$bwa_output/$sample.paired.bam"
    #sort the previous output
    samtools sort "$bwa_output/$sample.paired.bam" -o "$bwa_output/$sample.sorted.bam"
    #base recalibration trough the construction ofa data table 
    gatk BaseRecalibrator -I "$bwa_output/$sample.sorted.bam" -R ${ref_dir}/hg38.fasta --known-sites "$supporting_files_dir/Homo_sapiens_assembly38.dbsnp138.vcf" -O "$data_dir/$sample.data.table"
    #application of the previous generated data table
    gatk ApplyBQSR -I "$bwa_output/$sample.sorted.bam" -R ${ref_dir}/hg38.fasta --bqsr-recal-file "$data_dir/$sample.data.table" -O "$bwa_output/$sample.paired_sorted_bqsr.bam"
    #variants call with HaplotypeCaller
    gatk HaplotypeCaller -R "$ref_dir/hg38.fasta" -I "$bwa_output/$sample.paired_sorted_bqsr.bam" -O "$output_dir/$sample.vcf"

done
