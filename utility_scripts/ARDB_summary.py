__author__ = 'ksimmon'

import argparse
import sys

parser = argparse.ArgumentParser(description="the script take multiple jellyfish counts and compares the strains")
parser.add_argument("-i", "--input", help="cluster_output", type=argparse.FileType('r'), required=True)


args = parser.parse_args()


ARDB_file = args.input #"/data2/wash_u_sra/ENTEROCOCCUS/results/A_enterococcus_plasmid_cluster.txt"

results = {} #ardb code
current_strain = ""
ardb_code = ""
kmer_count = ""
ardb_info = ""

is_ardb_info = False
for line in ARDB_file:
    if "[ARDB INFORMATION]" in line:
        is_ardb_info = True
    if "[ARDB INFORMATION END]" in line:
        is_ardb_info = False

    if is_ardb_info and "[ARDB INFORMATION" not in line:

        if "STRAIN_ID:" in line or ":" not in line:
            current_strain = line.strip().split("\t")[-1].split("_")[0]

        elif "ARDB_CODE:" in line:
            ardb_code = line.strip().split("\t")[-1]

        elif "KMER_COUNT:" in line:
            kmer_count = line.strip().split("\t")[-1]

        elif "ARDB_INFO:" in line:
            ardb_info = line.strip().split("\t")[-1]
            if ardb_code in results:
                results[ardb_code].append((current_strain, kmer_count, ardb_info))
            else:
                results.update({ardb_code:[(current_strain,kmer_count,ardb_info)]})

for k, v in results.iteritems():
    sys.stdout.write("[{0}]\t{1}\n".format(k, v[0][2]))
    sys.stdout.write("strain\tkmer_count\n")
    for i in v:
        sys.stdout.write("{0}\t{1}\n".format(i[0],i[1]))
    print