import os, sys
import subprocess

base = "/Users/331-SimmonkLPTP/Documents/data/strain_typing/"

genome_file = "genomes_proks.txt"

for line in open(os.path.join(base,genome_file), "r"):
    if line[0] != "#":
        line = (line.strip().split("\t"))
        name = "_".join(line[0].split(" ")[0:2])
        strain = line[1].replace(" ", "_")
        clade = line[2].replace(" ", "_")
        refseq_FTP = line[-1]

        out = "|".join([name, strain, os.path.basename(refseq_FTP)])
        print(out)
        subprocess.check_call(["curl", "-o", os.path.join(base, out + ".fa.gz"),
                                "{0}/{1}_genomic.fna.gz".format(refseq_FTP, os.path.basename(refseq_FTP))])

