#!#!/home/ksimmon/anaconda/bin/python

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

def count_kmers(strain, path_to_kmers, queue):
    strain_objs[strain].kmer_count(path_to_kmers)



parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", help="List of jf databases seperated by commas (no_spaces)", type=str, )

args = parser.parse_args()
input_string = args.input

#####TESTING#########################################
_tmp = "E1_S20_clean_q25_31bp.jf,E2_S21_clean_q25_31bp.jf,E3_S22_clean_q25_31bp.jf,E4_S23_clean_q25_31bp.jf,EC_S24_clean_q25_31bp.jf"
base="/home/ksimmon/data/strain_typing/BactieraRun4_22_2014/enterococcus/jellyfish_31_clean/"
input_string = ""
for v in _tmp.split(","):
    input_string += "{0}{1},".format(base, v)
input_string = input_string[:-1]
####################################################

NUM_OF_STRAINS = len(input_string.split(","))

#strains = []
#file_paths = []

strain_objs = collections.OrderedDict()
for f in input_string.strip().split(","):
    #file_paths.append(f)
    #strains.append(os.path.basename(f).split("_")[0])
    stain_name = os.path.basename(f).split("_")[0]
    strain_objs.update({stain_name : jfobj.jf_object(stain_name, f)})

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
    pass#TODO

#create_dumped_file
try:
    subprocess.check_call(["jellyfish", "dump", "-tc",  "-L", "2", "-o", "/tmp/dump.fa", "/tmp/tmp.jf"] )
except:
    pass#TODO


count_table = []


count = 0

###start counting the shared content
q = Queue(maxsize=12)
for strain in strain_objs.iterkeys():
    p = Process(target=count_kmers, args=(strain, "/tmp/dump.fa", q))
    p.start()
    break

q.get()
p.join()


print strain_objs["E1"].coverage
print strain_objs["E1"].shared_count

count_table = []
for strain in strain_objs.itervalues():
    count_table.append(getattr(strain, "shared_count"))
    print getattr(strain, "shared_count")[1:10]


count_table = np.array(count_table)


print count_table[0].astype(bool)

print sum(count_table[0].astype(bool))
print sum(count_table[1].astype(bool))
print sum(count_table[0].astype(bool) == count_table[1].astype(bool))
print len(count_table[0])
print
print  float(sum(count_table[0].astype(bool) == count_table[1].astype(bool))) / float(len(count_table[0])) * 100



