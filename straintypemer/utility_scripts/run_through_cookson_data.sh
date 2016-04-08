#!/usr/bin/env bash

Acinetobacter_all="/home/ksimmon/storage/wash_u_sra/ACINETOBACTER/*jellyfish_31/*A*all*.jf"
Acinetobacter_plasmid="/home/ksimmon/storage/wash_u_sra/ACINETOBACTER/*jellyfish_31/*A*plasmid*.jf"
Acinetobacter_genome="~/storage/wash_u_sra/ACINETOBACTER/*jellyfish_31/*A*genome*.jf"

Acinetobacter_all_A="/home/ksimmon/storage/wash_u_sra/ACINETOBACTER/*jellyfish_31/*A_A*all*.jf"
Acinetobacter_plasmid_A="/home/ksimmon/storage/wash_u_sra/ACINETOBACTER/*jellyfish_31/*A_A*plasmid*.jf"
Acinetobacter_genome_A="/home/ksimmon/storage/wash_u_sra/ACINETOBACTER/*jellyfish_31/*A_A*genome*.jf"

Staphylococus_all="/home/ksimmon/storage/wash_u_sra/STAPHYLOCOCCUS/*jellyfish_31/*S*all*.jf"
Staphylococus_plasmid="/home/ksimmon/storage/wash_u_sra/STAPHYLOCOCCUS/*jellyfish_31/*S*plasmid*.jf"
Staphylococus_genome="/home/ksimmon/storage/wash_u_sra/STAPHYLOCOCCUS/*jellyfish_31/*S*genome*.jf"

Staphylococus_all_A="/home/ksimmon/storage/wash_u_sra/STAPHYLOCOCCUS/*jellyfish_31/*A_S*all*.jf"
Staphylococus_plasmid_A="/home/ksimmon/storage/wash_u_sra/STAPHYLOCOCCUS/*jellyfish_31/*A_S*plasmid*.jf"
Staphylococus_genome_A="/home/ksimmon/storage/wash_u_sra/STAPHYLOCOCCUS/*jellyfish_31/*A_S*genome*.jf"

Enterococcus_all="/home/ksimmon/storage/wash_u_sra/ENTEROCOCCUS/*jellyfish_31/*E*all*.jf"
Enterococcus_plasmid="/home/ksimmon/storage/wash_u_sra/ENTEROCOCCUS/*jellyfish_31/*E*plasmid*.jf"
Enterococcus_genome="/home/ksimmon/storage/wash_u_sra/ENTEROCOCCUS/*jellyfish_31/*E*genome*.jf"

Enterococcus_all_A="/home/ksimmon/storage/wash_u_sra/ENTEROCOCCUS/*jellyfish_31/*A_E*all*.jf"
Enterococcus_plasmid_A="/home/ksimmon/storage/wash_u_sra/ENTEROCOCCUS/*jellyfish_31/*A_E*plasmid*.jf"
Enterococcus_genome_A="/home/ksimmon/storage/wash_u_sra/ENTEROCOCCUS/*jellyfish_31/*A_E*genome*.jf"

python="/home/ksimmon/anaconda/bin/python"
StrainTypeMer="/home/ksimmon/PycharmProjects/strain_typing/Strain_typing/StrainTypeMer.py"
# 100,000 kmer cutoff
${python} ${StrainTypeMer} -c 100000 -t 8 --output_prefix ../temp_files/Acinetobacter_plasmid_A_100000 ${Acinetobacter_plasmid_A}
${python} ${StrainTypeMer} -c 100000 -t 8 --output_prefix ../temp_files/Acinetobacter_genome_A_100000 ${Acinetobacter_genome_A}
${python} ${StrainTypeMer} -c 100000 -t 8 --output_prefix ../temp_files/Acinetobacter_all_A_100000 ${Acinetobacter_all_A}

# 100,000 kmer cutoff
${python} ${StrainTypeMer} -t 8 --output_prefix ../temp_files/Acinetobacter_plasmid_A ${Acinetobacter_plasmid_A}
${python} ${StrainTypeMer} -t 8 --output_prefix ../temp_files/Acinetobacter_genome_A ${Acinetobacter_genome_A}
${python} ${StrainTypeMer} -t 8 --output_prefix ../temp_files/Acinetobacter_all_A ${Acinetobacter_all_A}

# 100,000 kmer cutoff
${python} ${StrainTypeMer} -c 100000 -t 8 --output_prefix ../temp_files/Acinetobacter_plasmid_100000 ${Acinetobacter_plasmid}
${python} ${StrainTypeMer} -c 100000 -t 8 --output_prefix ../temp_files/Acinetobacter_genome_100000 ${Acinetobacter_genome}
${python} ${StrainTypeMer} -c 100000 -t 8 --output_prefix ../temp_files/Acinetobacter_all_100000 ${Acinetobacter_all}

# 100,000 kmer cutoff
${python} ${StrainTypeMer} -t 8 --output_prefix ../temp_files/Acinetobacter_plasmid ${Acinetobacter_plasmid}
${python} ${StrainTypeMer} -t 8 --output_prefix ../temp_files/Acinetobacter_genome ${Acinetobacter_genome}
${python} ${StrainTypeMer} -t 8 --output_prefix ../temp_files/Acinetobacter_all ${Acinetobacter_all}

#########################################################################################################################

#cutoff
${python} ${StrainTypeMer} -c 100000 -t 8 --output_prefix ../temp_files/Staphylococus_plasmid_A_100000 ${Staphylococus_plasmid_A}
${python} ${StrainTypeMer} -c 100000 -t 8 --output_prefix ../temp_files/Staphylococus_genome_A_100000 ${Staphylococus_genome_A}
${python} ${StrainTypeMer} -c 100000 -t 8 --output_prefix ../temp_files/Staphylococus_all_A_100000 ${Staphylococus_all_A}

# 100,000 kmer cutoff
${python} ${StrainTypeMer} -t 8 --output_prefix ../temp_files/Staphylococus_plasmid_A ${Staphylococus_plasmid_A}
${python} ${StrainTypeMer} -t 8 --output_prefix ../temp_files/Staphylococus_genome_A ${Staphylococus_genome_A}
${python} ${StrainTypeMer} -t 8 --output_prefix ../temp_files/Staphylococus_all_A ${Staphylococus_all_A}

# 100,000 kmer cutoff
${python} ${StrainTypeMer} -c 100000 -t 8 --output_prefix ../temp_files/Staphylococus_plasmid_100000 ${Staphylococus_plasmid}
${python} ${StrainTypeMer} -c 100000 -t 8 --output_prefix ../temp_files/Staphylococus_genome_100000 ${Staphylococus_genome}
${python} ${StrainTypeMer} -c 100000 -t 8 --output_prefix ../temp_files/Staphylococus_all_100000 ${Staphylococus_all}

# 100,000 kmer cutoff
${python} ${StrainTypeMer} -t 8 --output_prefix ../temp_files/Staphylococus_plasmid ${Staphylococus_plasmid}
${python} ${StrainTypeMer} -t 8 --output_prefix ../temp_files/Staphylococus_genome ${Staphylococus_genome}
${python} ${StrainTypeMer} -t 8 --output_prefix ../temp_files/Staphylococus_all ${Staphylococus_all}

#########################################################################################################################

#cutoff
${python} ${StrainTypeMer} -c 100000 -t 8 --output_prefix ../temp_files/Enterococcus_plasmid_A_100000 ${Enterococcus_plasmid_A}
${python} ${StrainTypeMer} -c 100000 -t 8 --output_prefix ../temp_files/Enterococcus_genome_A_100000 ${Enterococcus_genome_A}
${python} ${StrainTypeMer} -c 100000 -t 8 --output_prefix ../temp_files/Enterococcus_all_A_100000 ${Enterococcus_all_A}

# 100,000 kmer cutoff
${python} ${StrainTypeMer} -t 8 --output_prefix ../temp_files/Enterococcus_plasmid_A ${Enterococcus_plasmid_A}
${python} ${StrainTypeMer} -t 8 --output_prefix ../temp_files/Enterococcus_genome_A ${Enterococcus_genome_A}
${python} ${StrainTypeMer} -t 8 --output_prefix ../temp_files/Enterococcus_all_A ${Enterococcus_all_A}

# 100,000 kmer cutoff
${python} ${StrainTypeMer} -c 100000 -t 8 --output_prefix ../temp_files/Enterococcus_plasmid_100000 ${Enterococcus_plasmid}
${python} ${StrainTypeMer} -c 100000 -t 8 --output_prefix ../temp_files/Enterococcus_genome_100000 ${Enterococcus_genome}
${python} ${StrainTypeMer} -c 100000 -t 8 --output_prefix ../temp_files/Enterococcus_all_100000 ${Enterococcus_all}

# 100,000 kmer cutoff
${python} ${StrainTypeMer} -t 8 --output_prefix ../temp_files/Enterococcus_plasmid ${Enterococcus_plasmid}
${python} ${StrainTypeMer} -t 8 --output_prefix ../temp_files/Enterococcus_genome ${Enterococcus_genome}
${python} ${StrainTypeMer} -t 8 --output_prefix ../temp_files/Enterococcus_all ${Enterococcus_all}

#No cleanup
