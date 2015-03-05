__author__ = 'keith simmon'

import numpy as np
import matplotlib.pyplot as plt


def parse_values(_file):
    _arr_count = []
    _arr_noKMERS = []
    with open(_file) as f:
        for line in f:
            line = line.strip().split(" ")
            _arr_count.append(float(line[0]))
            _arr_noKMERS.append(float(line[1]))


    return np.array(_arr_count), np.array(_arr_noKMERS)


"""
Exploring ways to qc NGS run to explore the approriate cutoff for the kmer count and determine need coverage.

Initially AC5 and AC4, which have been run 3 times each.

"""

from os import listdir
from os.path import isfile, join
base_dir = "/home/ksimmon/data/strain_typing/Acineto_combined/jellyfish_31/histo/"

files = listdir(base_dir)


#files=["A4_S4_31bp_31bp_histo.txt"]#, "AC9_2_S5_31bp_31bp_histo.txt", "A9_S9_31bp_31bp_histo.txt" ]



line = []
for f in files:
    if "6055" not in f:
        #if "A4" in f or "AC4" in f:
        print f
        x,y = parse_values(base_dir + f)
        #plt.bar(x,y, linewidth=1.0)

        n, bins, patches = plt.hist(x[1:], 50, normed=1, color="green", alpha=0.8)
        print np.mean(x), np.min(x),np.max(x)
        print np.mean(y)

#plt.bar(x,y )

        #plt.yscale("log")
        plt.ylabel("no of kmers")
        plt.xlabel("count")
        #plt.xlim(0,5000)
        #plt.ylim(0,5000)
        plt.show()


print 100 - 21 + 1
print 200 - 30 + 1
