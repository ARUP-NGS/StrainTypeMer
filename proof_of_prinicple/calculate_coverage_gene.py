__author__ = 'ksimmon'


from itertools import izip
import numpy as np
import math


file_2 = "A8_coverage.bed"
file_1 = "AC8_coverage.bed"


base = "/home/ksimmon/data/strain_typing/reference_strains/acinetobacter_genes/"



AC5_results = []
A52_results = []
gene = []

count = 1

for l1, l2 in izip(open(base+file_1), open(base+file_2)):
    A52_results.append(float(l1.strip().split("\t")[-1]))
    AC5_results.append(float(l2.strip().split("\t")[-1]))
    l = l2.split("\t")
    gene.append("{0}-{1}:{2}".format(count,l[1],l[2]))
    count+=1
    #print l

np_AC5 = np.array(AC5_results)
np_A52 = np.array(A52_results)

# print np.average(AC5_results), np.max(AC5_results), np.min(AC5_results)
# print np.average(A52_results), np.max(A52_results), np.min(A52_results)
# print np.sum(AC5_results)
# print np.sum(A52_results)

#print  (np.sum(AC5_results)) / np_AC5
#print  (np.sum(A52_results)) / np_A52

difference_dict = {}
values_dict = {}
print
for i in range(len(np_AC5)):
    r_A52 = 0
    r_AC5 = 0
    if np_A52[i] != 0:
        r_A52 = int(np.sum(A52_results)/ np_A52[i])
    if np_AC5[i] != 0:
        r_AC5 = int(np.sum(AC5_results)/ np_AC5[i])
    #print gene[i], r_AC5, r_A52
    values_dict.update({gene[i]:(r_AC5,r_A52)})
    difference_dict.update({gene[i]: abs(r_A52-r_AC5)})


for k in sorted(difference_dict, key=difference_dict.get, reverse=True):
    _tmp = k.split(":")
    length = int(_tmp[-1]) - int(_tmp[0].split("-")[-1])
    print "{0}\t{1}\t{2}\t{3}".format(difference_dict[k], k, length,  "\t".join(str(i) for i in values_dict[k]))





