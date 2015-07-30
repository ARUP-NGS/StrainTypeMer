__author__ = 'ksimmon'

import os, sys
base = os.path.abspath("/home/ksimmon/reference/mlst/acinetobacter/oxf")
from Bio import SeqIO
from Bio.Seq import Seq

reference = SeqIO.parse("/home/ksimmon/reference/mlst/acinetobacter/NC_021726.fa", "fasta").next()
wf = sys.stdout
wf = open("/home/ksimmon/reference/mlst/acinetobacter/NC_021726_Oxf.fa", "w")
#print len(reference)
pad = 200

def lower_padding(pad, sequence):
    return "{0}{1}{2}".format(sequence[0:200].lower(), sequence[200:-200], sequence[-200:].lower())

sequences = {}
for f in os.listdir(base):
    for s in SeqIO.parse(os.path.join(base,f), "fasta"):
        sequences.update({"Oxf_{0}".format(s.name.split("_")[0]) : s.seq   })
        break

for line in open("/home/ksimmon/reference/mlst/acinetobacter/oxf_acinetobacter_baumannii_NC_021726_sorted.bed"):
    ref, start, end, name, length, strand = line.strip().split("\t")
    if strand == "+":
        wf.write(">{0}\n{1}\n".format(name, lower_padding(pad,
                            str(reference.seq[int(start)-pad:int(end)+pad]))))
    else:
        wf.write(">{0}\n{1}\n".format(name, lower_padding(pad,
                            str(reference.seq[int(start)-pad:int(end)+pad].reverse_complement()))))
