from Bio import SeqIO
import os
import gzip
import subprocess
import sys

"""
this script is just going to calculate some basic info about the genomes at the jellyfish and fasta level

Genome
1: acc_number = genome size
2: replicons = size

unique kmers
total kmers
"""


def get_stats(jf_file):
    jellyfish_path = "/Users/331-SimmonkLPTP/bin/jellyfish-2.2.5/bin/jellyfish"
    op = subprocess.Popen([jellyfish_path, "stats", jf_file], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, err = op.communicate()
    if err != "":
        sys.stderr.write("Failed to retrieve db stats\n JELLYFISH ERROR: {0}\n".format(err))
        raise RuntimeError
    out = out.split("\n")
    unique_kmers = int(out[0].split(" ")[-1])
    distinct_kmers = int(out[1].split(" ")[-1])
    total_kmers = int(out[2].split(" ")[-1])
    max_count = int(out[3].split(" ")[-1])
    return unique_kmers, distinct_kmers, total_kmers, max_count


base = "/Users/331-SimmonkLPTP/Documents/data/strain_typing/acinetobacter_baumannii"
genome_path = os.path.join(base, "genomes")
counts_path = os.path.join(base, "counts_31")
for f in os.listdir(genome_path):
    if f[0] != ".":
        os.rename(os.path.join(genome_path, f), os.path.join(genome_path, f.replace(":","-")))
#
#         with gzip.open(os.path.join(genome_path, f)) as h:
#             for s in SeqIO.parse(h, "fasta"):
#                 print f, s.name, len(s)
#
# for f in os.listdir(counts_path):
#     print get_stats(os.path.join(counts_path, f))
