#!/usr/bin/env bash


cd ~/reference/strian_typing_resources/bacteria_plasmids/
files=*.fa


#efetch -id (basename $f .fa) -mode xml -db nucleotide |  xtract -element BinomialOrgName_genus BinomialOrgName_species >> plasmid_info_2.txt
for f in $files
    do
        echo "${f}"
        printf "${f%.fa}\t" >> ${f%.fa}_info.txt
        efetch -id "${f%.fa}" -mode xml -db nucleotide |  xtract -element BinomialOrgName_genus BinomialOrgName_species >> ${f%.fa}_info.txt &
    done
    