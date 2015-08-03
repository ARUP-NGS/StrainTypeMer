#!/usr/bin/env bash
###############################################################################################
# This bash script ties together reference assembly for MLST information
# @input fastq1
#         fastq2
#         reference_file
# 1) cutadapt - remove adapters
# 2) prinseq  - quality filtering
# 3) bwa      - mapped to reference (whole genome)
# 4) sort, index, pileup
# 5)
#############################################################################################

################   PATHS  ############################################################
#fastq_dir="/home/ksimmon/data/strain_typing/Acineto_combined/raw_fastqs/"
#mlst_dir="/home/ksimmon/data/strain_typing/Acineto_combined/MLST/"


reference_file="/home/ksimmon/reference/mlst/acinetobacter/NC_021726_Oxf.fa"

#### FILES
file_1_path=$1
file_2_path=$2
#reference_file=$3




file_1_name="${file_1_path##*/}"
file_2_name="${file_2_path##*/}"
sample_id=${file_1_name%_1.fastq.gz}
preprocessed_dir="${PWD}/${sample_id}_processed_pileup/"
log_file=${preprocessed_dir}${sample_id}.log

### check to make sure files and naming looks oK
if [ "${file_1_path}" == "${file_2_path}" ]; then
    echo "files have the same name"
    exit 1
fi

if [ "${sample_id}" != "${file_2_name%_2.fastq.gz}" ]; then
    echo "file names conflict"
    echo "${sample_id}"
    echo "${file_2_name%_2.fastq.gz}"
    exit 1
fi

if  [ ! -f ${file_1_path} ]; then
    echo "file not present ${file_1_path}"
    exit 1
fi

if  [ ! -f ${file_2_path} ]; then
    echo "file not present ${file_2_path}"
    exit 1
fi
#####################################################


#### APPLICATIONS #############################################################
cutadapt="/home/ksimmon/anaconda/bin/cutadapt"
prinseq="perl /home/ksimmon/bin/prinseq-lite-0.20.4/prinseq-lite.pl"
bwa="/home/ksimmon/bin/bwa-0.7.12/bwa"
samtools="/home/ksimmon/bin/samtools-1.2/samtools"
pileup_to_consensus="python /home/ksimmon/PycharmProjects/strain_typing/Strain_typing/utility_scripts/pileup_to_consensus.py"
bowtie2="/home/ksimmon/tools/bowtie2"
query_pubmlst="python /home/ksimmon/PycharmProjects/strain_typing/Strain_typing/utility_scripts/Query_pubmlst.py"


### VARIBLES ###########################################################
ADAPTER="AGATCGGAAGAGC"

##SETUP##
source /home/ksimmon/.bashrc

#cutadapt
mkdir ${preprocessed_dir}
touch ${log_file}
file_1_trimmed="${preprocessed_dir}${sample_id}_trimmed_1.fq"
file_2_trimmed="${preprocessed_dir}${sample_id}_trimmed_2.fq"


$cutadapt -a ${ADAPTER} -A ${ADAPTER} -o ${file_1_trimmed} -p ${file_2_trimmed} ${file_1_path} ${file_2_path} --minimum-length 1 >> ${log_file} 2>&1

######## PRINSEQ ###################
file_clean_base="${preprocessed_dir}${sample_id}_preprocessed"
file_failed_base="${preprocessed_dir}${sample_id}_failed"
file_clean_1="${preprocessed_dir}${sample_id}_preprocessed_1.fastq"
file_clean_2="${preprocessed_dir}${sample_id}_preprocessed_2.fastq"
file_clean_1_s="${preprocessed_dir}${sample_id}_preprocessed_1_singletons.fastq"
file_clean_2_s="${preprocessed_dir}${sample_id}_preprocessed_2_singletons.fastq"

prinseq_params="-out_format 3 -min_len 30 -ns_max_p 2 -lc_method dust -lc_threshold 10 -trim_tail_left 3 -trim_tail_right 3 -trim_qual_left 20 -trim_qual_right 20"
$prinseq -fastq ${file_1_trimmed} -fastq2 ${file_2_trimmed} -out_good ${file_clean_base} -out_bad ${file_failed_base} ${prinseq_params} >> ${log_file} 2>&1

##add singletons to files
file_1_processed="${preprocessed_dir}${sample_id}_prepro_1.fq"
file_2_processed="${preprocessed_dir}${sample_id}_prepro_2.fq"

###### BWA #####################
if  [ -s ${file_clean_1_s} ] && [ ${file_clean_2_s}  ]; then
    cat ${file_clean_1} ${file_clean_1_s} ${file_clean_2_s} > ${file_1_processed}
elif [ -s ${file_clean_1_s} ]; then
    cat ${file_clean_1} ${file_clean_1_s}  > ${file_1_processed}
elif [ -s ${file_clean_2_s} ]; then
    cat ${file_clean_1} ${file_clean_2_s} > ${file_1_processed}
else
    cat cat ${file_clean_1} > ${file_1_processed}

fi


cat ${file_clean_2} > ${file_2_processed}
## ALIGNMENT BWA ###
sam_file="${preprocessed_dir}${sample_id}.sam"
mapped_sam_file="${preprocessed_dir}${sample_id}_mapped.sam"
$bwa mem -t 8 ${reference_file} ${file_1_processed} ${file_2_processed} > ${sam_file} 2>> ${log_file}
#################################

#OR
######## BOWTIE #############
#$bowtie2 --local -p 8 -x ${reference_file} -1 ${file_clean_1} -2 ${file_clean_2} -U ${file_clean_1_s},${file_clean_2_s} -S ${sam_file}
#########bowtie2##################
#bioawk -c sam '!and($flag,4)' ${sam_file} > ${mapped_sam_file} #extract the mapped reads
bam_file="${sam_file%.sam}.bam"
$samtools view -b ${sam_file} > ${bam_file}

sorted_bam_file="${bam_file%.bam}_sorted"
$samtools sort ${bam_file} ${sorted_bam_file}

sorted_bam_file="${bam_file%.bam}_sorted.bam"
$samtools index "${sorted_bam_file}"

pile_up_file="${bam_file%.bam}.pileup"

$samtools mpileup -f ${reference_file} ${sorted_bam_file} > ${pile_up_file}


$pileup_to_consensus -r ${reference_file} -i ${pile_up_file} > ${preprocessed_dir}${sample_id}_consensus.fa

###get MLST STs and Alleles
$query_pubmlst  -i ${preprocessed_dir}${sample_id}_consensus.fa -d "${3}" > ${preprocessed_dir}${sample_id}_mlst.txt

echo "Done with ${sample_id}"
exit 0