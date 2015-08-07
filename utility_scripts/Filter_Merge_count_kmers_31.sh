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
kraken="/home/ksimmon/bin/kraken/kraken"

PATH="/home/ksimmon/bin/jellyfish-1.1.11/bin/:$PATH"

#### VARIBLES ###########################################################
ADAPTER="AGATCGGAAGAGC"

##SETUP##
#source /home/ksimmon/.bashrc

preprocessed_dir="${PWD}/${base_file_name}_processed_jellyfish_31/"
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
file_clean_1="${preprocessed_dir}${base_file_name}_preprocessed_1.fasta"
file_clean_2="${preprocessed_dir}${base_file_name}_preprocessed_2.fasta"
file_clean_1_s="${preprocessed_dir}${base_file_name}_preprocessed_1_singletons.fasta"
file_clean_2_s="${preprocessed_dir}${base_file_name}_preprocessed_2_singletons.fasta"

file_concat="${preprocessed_dir}${base_file_name}_preprocessed_1_2.fasta"
file_concat_all="${preprocessed_dir}${base_file_name}_preprocessed_all.fasta"

echo "quality filtering reads"
prinseq_params="-out_format 1 -line_width 0 -derep 1 -min_len 30 -ns_max_p 2 -lc_method dust -lc_threshold 10 -trim_tail_left 3 -trim_tail_right 3 -trim_qual_left 20 -trim_qual_right 20"
$prinseq -fastq ${file_1_trimmed} -fastq2 ${file_2_trimmed} -out_good ${file_clean_base} -out_bad ${file_failed_base} ${prinseq_params} >> ${log_file} 2>&1
#
###add singletons to files
#concatenate duplicates

$concat_fastq -1 ${file_clean_1} -2 ${file_clean_2} -o ${file_concat} -f fa

echo "combining reads"
cat ${file_concat} ${file_clean_1_s} ${file_clean_2_s} > ${file_concat_all}


##ADD KRAKEN PREPROCESSING
kraken_db=" --db /home/ksimmon/reference/strian_typing_resources/kraken_strain_typing/ "

kraken_options=" --threads 8 --fasta-input "

echo "identifying plasmid and genomic sequences "
$kraken ${kraken_db} ${kraken_options} ${file_concat_all} > "${preprocessed_dir}${base_file_name}_kraken_raw.txt"

plasmid_fasta="${preprocessed_dir}${base_file_name}_plasmid_sequence.fa"
genomic_fasta="${preprocessed_dir}${base_file_name}_genome_sequence.fa"

echo "sequestering reads"
python "/home/ksimmon/PycharmProjects/strain_typing/Strain_typing/utility_scripts/extract_reads.py" \
--kraken_input "${preprocessed_dir}${base_file_name}_kraken_raw.txt" \
--fasta_file "${file_concat_all}" \
--taxonomy_file "/home/ksimmon/reference/strian_typing_resources/kraken_strain_typing/taxonomy/nodes.dmp" \
--parent_taxid ${4} > "${genomic_fasta}" &

python "/home/ksimmon/PycharmProjects/strain_typing/Strain_typing/utility_scripts/extract_reads.py" \
--kraken_input "${preprocessed_dir}${base_file_name}_kraken_raw.txt" \
--fasta_file "${file_concat_all}" \
--taxonomy_file "/home/ksimmon/reference/strian_typing_resources/kraken_strain_typing/taxonomy/nodes.dmp" \
--parent_taxid 6 > "${plasmid_fasta}" &


wait
echo "counting kmers"
jellyfish_out="${preprocessed_dir}${base_file_name}_jellyfish_${3}.jf"

$jellyfish count -m ${3} -L 2 -t 10 -o ${jellyfish_out} -s 4G -C ${file_concat_all}

jellyfish_plasmid_out="${preprocessed_dir}${base_file_name}_plasmid_jellyfish_${3}.jf"
$jellyfish count -m ${3} -L 2 -t 10 -o ${jellyfish_plasmid_out} -s 4G -C ${plasmid_fasta}

jellyfish_genome_out="${preprocessed_dir}${base_file_name}_genome_jellyfish_${3}.jf"
$jellyfish count -m ${3} -L 2 -t 10 -o ${jellyfish_genome_out} -s 4G -C ${genomic_fasta}

echo "script finished"