#!/home/ksimmon/anaconda/bin/python

__author__ = 'ksimmon'


"""

This script will take as input jf databases from strains and perform a comparison

The input of strains should be provided should be provided as a comma seperated list.
assume sample name is before underscore

"""


import jellyfish
import numpy as np
import argparse
import os.path
import sys
import subprocess
import collections
import jf_object as jfobj
from multiprocessing import Process, Queue

def count_kmers(strain_name, path_to_kmers, queue):
    queue.put(strain_objs[strain_name].kmer_count(strain_name, path_to_kmers))
    #print _a[1:10]


parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", help="List of jf databases seperated by commas (no_spaces)", type=str,)

args = parser.parse_args()
input_string = args.input

#####TESTING#########################################
_tmp = "E1_S20_clean_q25_31bp.jf,E2_S21_clean_q25_31bp.jf"#,E3_S22_clean_q25_31bp.jf,E4_S23_clean_q25_31bp.jf,EC_S24_clean_q25_31bp.jf"
base="/home/ksimmon/data/strain_typing/BactieraRun4_22_2014/enterococcus/jellyfish_31_clean/"
#_tmp = "S1_S12_clean_q25_31bp.jf,S2_S13_clean_q25_31bp.jf,S3_S14_clean_q25_31bp.jf,S4_S15_clean_q25_31bp.jf,S5_S16_clean_q25_31bp.jf,S6_S17_clean_q25_31bp.jf,S7_S18_clean_q25_31bp.jf,SC_S19_clean_q25_31bp.jf"
#base = "/home/ksimmon/data/strain_typing/BactieraRun4_22_2014/staph/jellyfish_31_clean/"
input_string = ""
for v in _tmp.split(","):
    input_string += "{0}{1},".format(base, v)
input_string = input_string[:-1]
####################################################

NUM_OF_STRAINS = len(input_string.split(","))

#strains = []
file_paths = []

strain_objs = collections.OrderedDict()
for f in input_string.strip().split(","):
    file_paths.append(f)
    #strains.append(os.path.basename(f).split("_")[0])
    strain_name = os.path.basename(f).split("_")[0]
    strain_objs.update({strain_name : jfobj.jf_object(strain_name, f)})

    if os.path.exists(f) == False:
        sys.stderr.write("file does not exist: {0}\n".format(f))
        sys.exit(1)


if len(strain_objs) != NUM_OF_STRAINS:
    sys.stderr.write("strain names are not unique:\n{0}\n".format(", ".join(strain_objs.keys())))
    sys.exit(1)

#create merged file
try:
    subprocess.check_call(["jellyfish", "merge", "-o", "/tmp/tmp.jf"] + file_paths )
except:
    pass


#create_dumped_file
try:
    subprocess.check_call(["jellyfish", "dump", "-tc",  "-L", "2", "-o", "/tmp/dump.fa", "/tmp/tmp.jf"] )
except:
    pass#TODO


count_table = []


count = 0

###start counting the shared content
# q = Queue()
#
# for k in strain_objs.iterkeys():  #print k
#     p = Process(target=count_kmers, args=(k, "/tmp/dump.fa", q))
#     p.start()
# p.join()
#
#
# while not q.empty():
#     strain, counts = q.get()
#     setattr(strain_objs[strain], "shared_count", counts)
count_table = []
count = 0
with open("/tmp/dump.fa") as f:
    for l in f:
        _arr = []
        attach = _arr.append
        count += 1
        kmer = l.strip().split("\t")[0]
        mer = jellyfish.MerDNA(kmer)
        mer.canonicalize()
        for v in strain_objs.itervalues():
            #answers.append(str(query_db[mer]))
            #print
            attach(int(v.qf[mer]))
        if  np.array(_arr).min() == 0:
            if np.array(_arr).max() >4:
                print kmer, _arr


        count_table.append(_arr)


        if count > 10000:
            break

lowest_coverage = 10000

sys.stdout.write(" COVERAGE INFORMATION ".center(80, "-") + "\n")
sys.stdout.write("calculations only include kmers counted > 3 times\n")

cutoff = []
for k, v in strain_objs.iteritems():
    sys.stdout.write("Strain: {:s}\t Coverage Estimate: {:.1f}\n".format(k, v.coverage))
    cutoff.append(int(v.coverage * .20))
    v.set_cutoff(int(v.coverage * .20))
    if v.coverage < lowest_coverage:
        strain, lowest_coverage = k, v.coverage
sys.stdout.write("".center(80, "-") + "\n\n")

filter_kmer_count = int(lowest_coverage * .25)

count_table = np.array(count_table)
sys.stdout.write("{:s} has the lowest coverage, estimated at: {:.1f}\n".format(strain,lowest_coverage))
sys.stdout.write("Filtering out kmers with counts lower than {:.0f}\n\n".format(filter_kmer_count))



sys.stdout.write(" STRAIN STATS ".center(80, "-") + "\n")
for k, v in strain_objs.iteritems():
    sys.stdout.write("Strain: {:s}\n".format(k))
    sys.stdout.write("\tInferred genome size: {:.0f}  [filtering kmers counted fewer than {:.0f} times]\n".format(
                      v.estimate_genome_size(v.kmer_cutoff), v.kmer_cutoff ))
sys.stdout.write("".center(80, "-") + "\n\n")

cutoff = np.array(cutoff)
##FILTER COUNT TABLE
count_table_filtered = []
total_kmers = 0.0
kept_kmers = 0.0
for kmer_count in count_table:
    total_kmers += 1
    if (np.array(kmer_count) - cutoff).max() >= 0:
        count_table_filtered.append(kmer_count)
        kept_kmers += 1


sys.stdout.write(" FILTER INFO ".center(80, "-") + "\n")
sys.stdout.write("{:.0f}\tkmers counted (n > 1) in all strains\n".format(total_kmers))
sys.stdout.write("{:.0f} ({:.1f}%)\tkmers remaining after filter in all strains\n".format(kept_kmers, float(kept_kmers)/total_kmers * 100.0))
sys.stdout.write("".center(80, "-") + "\n\n")

count_table_filtered = np.array(count_table_filtered).T

#print count_table_filtered[4]

#print  float(sum(count_table_filtered[0].astype(bool) == count_table_filtered[1].astype(bool))) / float(len(count_table_filtered[0])) * 100


strain_keys = strain_objs.keys()

#DETERMINE RELATIONSHIPS


similarity_dict = collections.OrderedDict()
#_str = "\t" + "\t".join(strain_keys) + "\n"
for i in range(len(count_table_filtered) ):
    #_str += strain_keys[i] + "\t"
    similarity_dict.update({strain_keys[i] : collections.OrderedDict()})
    for j in range( i + 1, len(count_table_filtered)):
        similarity_dict[strain_keys[i]].update({strain_keys[j] :
            float(sum(count_table_filtered[i].astype(bool) == count_table_filtered[j].astype(bool))) /
            float(len(count_table_filtered[0])) * 100 })
    #_str += "\n"

_DEL = ","

_str = _DEL + _DEL.join(strain_keys) + "\n"
for i in range(len(strain_keys)):
    _str += strain_keys[i] + _DEL
    for j in range(len(strain_keys)):
        if i == j:
            _str +=  "100" + _DEL
        elif i < j:
            _str += "{:.1f}{:s}".format(similarity_dict[strain_keys[i]][strain_keys[j]],_DEL)
        elif i > j:
            _str += "{:.1f}{:s}".format(similarity_dict[strain_keys[j]][strain_keys[i]],_DEL)

    _str  = _str[:-1] + "\n"

print _str