#!/usr/bin/env bash

########################################################################################################################
#
# This script takes paired fastq files as input runs:
#
# cutadapt - for trimming adapters
# prinseq - to remove and trim bad quality sequence
# merge_pair script
# jellyfish
#
# $1 fastq_1.gz
# $2 fastq_2.gz
# $3 kmer_length (for jellyfish) [default = 31]
source /home/ksimmon/.bashrc

fastq_1=${1##*/} #get the filename without dir
fastq_2=${2##*/}

#echo ${fastq_1%_*}
#echo ${fastq_2%_*}

if [ ${fastq_1%_*} != ${fastq_2%_*} ]
then
    >&2 echo "base file names not equal"
    >&2 echo "make sure the files are properly mated"
    exit 1;
fi

base_file_name=${fastq_1%_*}


kmer_length=$3

cutadapt="/home/ksimmon/anaconda/bin/cutadapt"
prinseq="perl /home/ksimmon/bin/prinseq-lite-0.20.4/prinseq-lite.pl"
jellyfish="/home/ksimmon/bin/jellyfish-2.2.0/bin/jellyfish"
concat_fastq="python /home/ksimmon/bin/Taxonomer/utilities/concatenate_fastqs_2.0.py"



#### VARIBLES ###########################################################
ADAPTER="AGATCGGAAGAGC"

##SETUP##
#source /home/ksimmon/.bashrc

preprocessed_dir="${PWD}/${base_file_name}_preprocessed/"
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
$cutadapt -a ${ADAPTER} -A ${ADAPTER} -o ${file_1_trimmed} -p ${file_2_trimmed} ${1} ${2} --minimum-length 1 >> ${log_file} 2>&1
#
######### PRINSEQ ###################
file_clean_base="${preprocessed_dir}${base_file_name}_preprocessed"
file_failed_base="${preprocessed_dir}${base_file_name}_failed"
file_clean_1="${preprocessed_dir}${base_file_name}_preprocessed_1.fasta"
file_clean_2="${preprocessed_dir}${base_file_name}_preprocessed_2.fasta"
file_clean_1_s="${preprocessed_dir}${base_file_name}_preprocessed_1_singletons.fasta"
file_clean_2_s="${preprocessed_dir}${base_file_name}_preprocessed_2_singletons.fasta"

file_concat="${preprocessed_dir}${base_file_name}_preprocessed_1_2.fasta"
file_concat_all="${preprocessed_dir}${base_file_name}_preprocessed_all.fasta"
#
prinseq_params="-out_format 1 -line_width 0 -min_len 30 -ns_max_p 2 -lc_method dust -lc_threshold 10 -trim_tail_left 3 -trim_tail_right 3 -trim_qual_left 20 -trim_qual_right 20"
$prinseq -fastq ${file_1_trimmed} -fastq2 ${file_2_trimmed} -out_good ${file_clean_base} -out_bad ${file_failed_base} ${prinseq_params} >> ${log_file} 2>&1
#
###add singletons to files
#concatenate duplicates
$concat_fastq -1 ${file_clean_1} -2 ${file_clean_2} -o ${file_concat} -f fa

cat ${file_concat} ${file_clean_1_s} ${file_clean_2_s} > ${file_concat_all}

jellyfish_out="${preprocessed_dir}${base_file_name}_jellyfish_${3}.jf"
$jellyfish count -m $3 -t 10 -o ${jellyfish_out} -s 4G -C ${file_concat_all}