__author__ = 'ksimmon'


import os
import os.path
import sys

base="/home/ksimmon/data/strain_typing/wash_u_sra/"

fastqs =  [i for i in os.listdir(base) if i.endswith("fastq.gz")]
print fastqs

file_handles = {}
for i in open(base + "SraRunTable.txt"):
    i = i.split("\t")
    #print i
    fastq_base, sample_name, species = i[5], i[1], i[4]
    #print fastq_base
    file_handles.update({ fastq_base  : "{0}_{1}_{2}".format(sample_name, species.replace(" ","_"), fastq_base) })


for i in fastqs:
    SRR = i.split("_")[0]
    postfix = "_" + i.split("_")[1]


    print file_handles[SRR] + postfix
    os.rename(base + i, base + file_handles[SRR] + postfix)