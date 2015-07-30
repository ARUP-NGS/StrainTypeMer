__author__ = 'ksimmon'

_tmp = "E1_S20_60bp.fa"
#,E2_S21_60bp.jf"#,E3_S22_60bp.jf,E4_S23_60bp.jf,EC_S24_60bp.jf"
#base="/home/ksimmon/data/strain_typing/BactieraRun4_22_2014/enterococcus/jellyfish_31_clean/"
base="/home/ksimmon/data/strain_typing/BactieraRun4_22_2014/enterococcus/jellyfish_60_clean/"

from Bio import SeqIO
import numpy as np


count_dict = {}
for s in SeqIO.parse(base+_tmp, "fasta"):
    c = int(s.name)
    v = float((s.seq.count("A")) + (s.seq.count("T"))) / 60.0 * 100
    #+ s.seq.count("T")) / 60
    if c in count_dict:
        count_dict[c][0] += v
        count_dict[c][1] += 100.0
    else:
        count_dict.update({c:[v,100.0]})




for k in sorted(count_dict):
    print "{:.0f}\t{:.1f}%".format(k, float(count_dict[k][0]) / count_dict[k][1])