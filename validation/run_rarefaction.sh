#!/usr/bin/env bash

base_path="/Users/331-SimmonkLPTP/Box Sync/ARUP/ACI_fastqs/clean_files/"
files="A10A_trimmed.fastq.gz
A10B_trimmed.fastq.gz
A11A_trimmed.fastq.gz
A11B_trimmed.fastq.gz
A12A_trimmed.fastq.gz
A12B_trimmed.fastq.gz
A13A_trimmed.fastq.gz
A13B_trimmed.fastq.gz
A14A_trimmed.fastq.gz
A14B_trimmed.fastq.gz
A15A_trimmed.fastq.gz
A15B_trimmed.fastq.gz
A16A_trimmed.fastq.gz
A16B_trimmed.fastq.gz
A1A_trimmed.fastq.gz
A1B_trimmed.fastq.gz
A2A_trimmed.fastq.gz
A2B_trimmed.fastq.gz
A3A_trimmed.fastq.gz
A3B_trimmed.fastq.gz
A5A_trimmed.fastq.gz
A5B_trimmed.fastq.gz
A6A_trimmed.fastq.gz
A6B_trimmed.fastq.gz
A7A_trimmed.fastq.gz
A7B_trimmed.fastq.gz
A8A_trimmed.fastq.gz
A8B_trimmed.fastq.gz
A9A_trimmed.fastq.gz
A9B_trimmed.fastq.gz"

for f in $files
    do
    python ${HOME}/git_repos/StrainTypeMer/validation/rarefaction_measure -i "${base_path}${f}" -o "${base_path}${f%%_trimmed.fastq.gz}_rare_counts.txt" &
done
wait