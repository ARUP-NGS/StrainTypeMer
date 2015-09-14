#!/usr/bin/env bash

########################################################################################################################
#
# This script takes paired fastq files as input runs:
#
# cutadapt - for trimming adapters
# prinseq - to remove and trim bad quality sequence
#
#
# $1 fastq_1.gz
# $2 fastq_2.gz
# $3 base_file_name

source /home/ksimmon/.bashrc

fastq_1=$(basename $1) #get the filename without dir
fastq_2=$(basename $2)



base_file_name="${3}"
echo "Prefix: ${base_file_name}"
echo "file 1: ${fastq_1}"
echo "file 2: ${fastq_2}"
reference="${4}"
echo "reference_file: ${reference}"

cutadapt="/home/ksimmon/anaconda/bin/cutadapt"
prinseq="perl /home/ksimmon/bin/prinseq-lite-0.20.4/prinseq-lite.pl"
bwa="/home/ksimmon/bin/bwa-0.7.12/bwa"
samtools="/home/ksimmon/bin/samtools-1.2/samtools"
bcftools="/home/ksimmon/bin/bcftools/bcftools"


#reference="/home/ksimmon/reference/strian_typing_resources/acinetobacter_genomes/CU459141.fa"


#jellyfish="/home/ksimmon/bin/jellyfish-2.2.0/bin/jellyfish"
#concat_fastq="python /home/ksimmon/bin/Taxonomer/utilities/concatenate_fastqs_2.0.py"
#kraken="/home/ksimmon/bin/kraken/kraken"
#PATH="/home/ksimmon/bin/jellyfish-1.1.11/bin/:$PATH"

#### VARIBLES ###########################################################
ADAPTER="AGATCGGAAGAGC"

##SETUP##
#source /home/ksimmon/.bashrc

preprocessed_dir="${PWD}/${base_file_name}_cookson_procedure/"
##cutadapt
mkdir ${preprocessed_dir}
log_file="${preprocessed_dir}${base_file_name}.log"
#
#
#
file_1_trimmed="${preprocessed_dir}${base_file_name}_trimmed_1.fq"
file_2_trimmed="${preprocessed_dir}${base_file_name}_trimmed_2.fq"
#
#
echo "trimming adapters"
$cutadapt -a ${ADAPTER} -A ${ADAPTER} -o ${file_1_trimmed} -p ${file_2_trimmed} ${1} ${2} --minimum-length 1 >> ${log_file} 2>&1
#
######### PRINSEQ ###################
file_clean_base="${preprocessed_dir}${base_file_name}_preprocessed"
file_failed_base="${preprocessed_dir}${base_file_name}_failed"
file_clean_1="${preprocessed_dir}${base_file_name}_preprocessed_1.fastq"
file_clean_2="${preprocessed_dir}${base_file_name}_preprocessed_2.fastq"
file_clean_1_s="${preprocessed_dir}${base_file_name}_preprocessed_1_singletons.fastq"
file_clean_2_s="${preprocessed_dir}${base_file_name}_preprocessed_2_singletons.fastq"


echo "quality filtering reads"
prinseq_params=" -out_format 3 -line_width 0 -derep 1 -min_len 30 -ns_max_p 2 "
$prinseq  -fastq ${file_1_trimmed} -fastq2 ${file_2_trimmed} -out_good ${file_clean_base} -out_bad ${file_failed_base} ${prinseq_params} >> ${log_file} 2>&1
#
sam_file="${preprocessed_dir}${base_file_name}.sam"

$bwa mem -t 4 -T 10 "${reference}" "${file_clean_1}" "${file_clean_2}" > "${sam_file}"


bam_file="${preprocessed_dir}${base_file_name}.bam"
$samtools view -b ${sam_file} > ${bam_file}

sorted_bam_file="${bam_file%.bam}_sorted"
$samtools sort ${bam_file} ${sorted_bam_file}

sorted_bam_file="${bam_file%.bam}_sorted.bam"
$samtools index "${sorted_bam_file}"

pile_up_file="${bam_file%.bam}.pileup"
$samtools mpileup -f "${reference}" "${sorted_bam_file}" > "${pile_up_file}"

vcf_file="${bam_file%.bam}.vcf"
$samtools mpileup -t "DP,DPR,DV,DP4,INFO/DPR,SP" -uf "${reference}" "${sorted_bam_file}" | $bcftools call -m -O v > "${vcf_file}"

echo "done ${base_file_name}"
