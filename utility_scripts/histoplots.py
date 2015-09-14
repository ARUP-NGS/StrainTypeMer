__author__ = 'ksimmon'

import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import argparse

parser = argparse.ArgumentParser(description="the script take multiple jellyfish counts and compares the strains")
parser.add_argument('jf_histo_files', nargs='+', help='jellyfish files for each strain', type=argparse.FileType("r") )

args = parser.parse_args()
jf_files = args.jf_histo_files

complete_results = []
for f in jf_files:
    _d = {i+1: 0  for i in range(300)}
    for line in open(f):
        line = line.strip().split(" ")
        freq, count = [int(i) for i in line]
        if freq in _d:
            _d[freq]+=count
        else:
            _d[300]+=count
    complete_results.append( (f.split("/")[1].split("_")[0],_d,) )

for name,_d in complete_results:
    N = len(_d)
    freq = _d.keys()
    count =   _d.values()
    ind = np.arange(N)  # the x locations for the groups
    width = 1      # the width of the bars
    fig, ax = plt.subplots(figsize=(20,4))
    #plt.figure(figsize=(3,4))
    rects1 = ax.bar(freq, count, width, color='r')
    ax.set_title(name + ": histogram of kmers observed at different counts")
    ax.set_xlabel("Frequency at which a kmer occurs")
    ax.set_ylabel("Distinct kmers")
    plt.xlim(0,150)
    plt.ylim(0,350000)
    plt.savefig(name + "_histo.png",)
#plt.show()