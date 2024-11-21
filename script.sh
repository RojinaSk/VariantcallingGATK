#!/bin/bash

# Script to call germline variants in a human WGS paired end reads 2 X 100bp
# Following GATK4 best practices workflow - https://gatk.broadinstitute.org/hc/en-us/articles/360035535932-Germline-short-variant-discovery-SNPs-Indels-
# This script is for demonstration purposes only

if false
then
# # download data from the 1000 genomes project
wget -P /users/5/sapko041/Ind_projects/variant_calling_GATK4/reads ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/phase3/data/HG00096/sequence_read/SRR062634_1.filt.fastq.gz

wget -P /users/5/sapko041/Ind_projects/variant_calling_GATK4/reads ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/phase3/data/HG00096/sequence_read/SRR062634_2.filt.fastq.gz


echo "Run Prep files..."

################################################### Prep files (TO BE GENERATED ONLY ONCE) ##########################################################



# download reference files and unzip it
 wget -P ~/users/5/sapko041/Ind_projects/variant_calling_GATK4/supporting_files/hg38/ https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
 gunzip ~/users/5/sapko041/Ind_projects/variant_calling_GATK4/supporting_files/hg38/hg38.fa.gz

# index ref - .fai file before running haplotype caller

samtools faidx ~/users/5/sapko041/Ind_projects/variant_calling_GATK4/supporting_files/hg38/hg38.fa


# ref dict - .dict file before running haplotype caller
gatk CreateSequenceDictionary R=~/users/5/sapko041/Ind_projects/variant_calling_GATK4/supporting_files/hg38/hg38.fa O=~/users/5/sapko041/Ind_projects/variant_calling_GATK4/supporting_files/hg38/hg38.dict


# download known sites files for BQSR from GATK resource bundle
wget -P ~/users/5/sapko041/Ind_projects/variant_calling_GATK4/supporting_files/hg38/ https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf
wget -P ~/users/5/sapko041/Ind_projects/variant_calling_GATK4/supporting_files/hg38/ https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf.idx



###################################################### VARIANT CALLING STEPS ####################################################################
fi

# directories
#saving paths for easier downstream analysis
ref="/users/5/sapko041/Ind_projects/variant_calling_GATK4/supporting_files/hg38/hg38.fa"
known_sites="/users/5/sapko041/Ind_projects/variant_calling_GATK4/supporting_files/hg38/Homo_sapiens_assembly38.dbsnp138.vcf"
aligned_reads="/users/5/sapko041/Ind_projects/variant_calling_GATK4/aligned_reads"
reads="/users/5/sapko041/Ind_projects/variant_calling_GATK4/reads"
results="/users/5/sapko041/Ind_projects/variant_calling_GATK4/results"
data="/users/5/sapko041/Ind_projects/variant_calling_GATK4/data"




# -------------------
# STEP 1: QC - Run fastqc 
# -------------------since paired end reads, we have two files

#echo "STEP 1: QC - Run fastqc"

fastqc ${reads}/SRR062634_1.filt.fastq.gz -o ${reads}/ #contains forward pairs
fastqc ${reads}/SRR062634_2.filt.fastq.gz -o ${reads}/ #reverse reads

# # No trimming required, quality looks okay.


# # --------------------------------------
# # STEP 2: Map to reference using BWA-MEM
# # --------------------------------------

  echo "STEP 2: Map to reference using BWA-MEM"

# # # BWA index reference(generate index for the reference)
 bwa index ${ref}


# # # BWA alignment using bwa mem
#read group info in " ". tPL(platform information), tSM(sample information) ${reference file} in reads folder > output the aligned reads in aligned reads folder
 bwa mem -t 4 -R "@RG\tID:SRR062634\tPL:ILLUMINA\tSM:SRR062634" ${ref} ${reads}/SRR062634_1.filt.fastq.gz ${reads}/SRR062634_2.filt.fastq.gz > ${aligned_reads}/SRR062634.paired.sam




# # -----------------------------------------
# # STEP 3: Mark Duplicates and Sort - GATK4
# # -----------------------------------------
#performs flagging and sorting the sam files
echo "STEP 3: Mark Duplicates and Sort - GATK4"
#provide the aligned reads file and output to aligned reads as sorted file
gatk MarkDuplicatesSpark -I ${aligned_reads}/SRR062634.paired.sam -O ${aligned_reads}/SRR062634_sorted_dedup_reads.bam



# # ----------------------------------
# # STEP 4: Base quality recalibration
# # ----------------------------------

#build model usimng GATK base recalibrator function
echo "STEP 4: Base quality recalibration"

# # 1. build the model(-I inpur, -R reference file -known sites -o output a table in data folder)
gatk BaseRecalibrator -I ${aligned_reads}/SRR062634_sorted_dedup_reads.bam -R ${ref} --known-sites ${known_sites} -O ${data}/recal_data.table


# # 2. Apply the model to adjust the base quality scores(ApplyBQSR)
gatk ApplyBQSR -I ${aligned_reads}/SRR062634_sorted_dedup_reads.bam -R ${ref} --bqsr-recal-file {$data}/recal_data.table -O ${aligned_reads}/SRR062634_sorted_dedup_bqsr_reads.bam 



# # -----------------------------------------------
# # STEP 5: Collect Alignment & Insert Size Metrics(QC)
# # -----------------------------------------------
#the I and INPUT diff is based on the clarity and brevity and also based on the tools used in GATK

echo "STEP 5: Collect Alignment & Insert Size Metrics"

gatk CollectAlignmentSummaryMetrics R=${ref} I=${aligned_reads}/SRR062634_sorted_dedup_bqsr_reads.bam O=${aligned_reads}/alignment_metrics.txt
gatk CollectInsertSizeMetrics INPUT=${aligned_reads}/SRR062634_sorted_dedup_bqsr_reads.bam OUTPUT=${aligned_reads}/insert_size_metrics.txt HISTOGRAM_FILE=${aligned_reads}/insert_size_histogram.pdf



# # ----------------------------------------------
# # STEP 6: Call Variants - gatk haplotype caller
# # ----------------------------------------------
#calling variants on analysis ready reads
echo "STEP 6: Call Variants - gatk haplotype caller"

gatk HaplotypeCaller -R ${ref} -I ${aligned_reads}/SRR062634_sorted_dedup_bqsr_reads.bam -O ${results}/raw_variants.vcf



# # extract SNPs & INDELS
#-R reference file, -V variants file , -O output, --select type (type of variant)
gatk SelectVariants -R ${ref} -V ${results}/raw_variants.vcf --select-type SNP -O ${results}/raw_snps.vcf
gatk SelectVariants -R ${ref} -V ${results}/raw_variants.vcf --select-type INDEL -O ${results}/raw_indels.vcf






