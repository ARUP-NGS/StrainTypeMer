__author__ = 'ksimmon'

def parse_pileup(coverage=0,bases="", quals="", line=""):
    quals = [ord(i)-33 for i in quals]
    if len(quals) != coverage:
        sys.stderr.write("The reads number and number of quality scores does not agree")
        sys.exit(1)
    bases = bases.replace("$", "").replace(",",".") #remove line starts and , to .
    enum = enumerate(bases)
    _bases = []
    for i,v in enum:
        if v is '^':
            enum.next()
            i,v = enum.next()
        if v in ["-", "+"]:
            _s = v
            _d = ""
            indel = True
            while indel:
                i,v = enum.next()
                while v.isdigit():
                    _d += v
                    i,v = enum.next()
                _s += _d
                _d = int(_d)
                for j in range(_d):
                    _s += v
                    i,v = enum.next()
                indel = False
            _bases.append(_s)
            _bases.append(v)
        else:
            _bases.append(v)
    print _bases
    print quals
    return _bases, quals

########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################


from Bio import SeqIO
import sys
import argparse
from collections import Counter

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", help="mpile file for samtools", type=str,)# required=True)
parser.add_argument("-r", "--reference_file", help="padded reference files Ns will be stripped",
                    type=str,)# required=True)
args = parser.parse_args()
pu_file = args.input
reference_file = args.reference_file

if pu_file == None and reference_file == None:
    pass
    ##TESTING
    pu_file = "/home/ksimmon/data/strain_typing/MiSeq_run_July_2015/renamed/148-7_3-Acinetobacter-gp3_processed_pileup/148-7_3-Acinetobacter-gp3.pileup"
    reference_file = "/home/ksimmon/reference/mlst/acinetobacter/NC_021726_Oxf.fa"

reference = {}
reference_consensus = {}
for s in SeqIO.parse(reference_file, "fasta"):
    #print s.name
    ind5 = -1
    ind3 = -1
    for i,char in enumerate(str(s.seq)):
        if char.isupper():
            if ind5 is -1:
                ind5 = i
            ind3 = i + 1
    reference.update({s.name:(ind5,ind3)})
    reference_consensus.update({s.name: ""})
#print reference

for line in open(pu_file, "r"):
    l = line.strip().split()
    ref_name, pos, ref, num_reads = l[0], int(l[1]), l[2], int(l[3])

    if pos >= reference[ref_name][0] - 10 and pos <= reference[ref_name][1] + 10:
        #print l
        try:
            bases, quals = parse_pileup(coverage=num_reads, bases=l[4], quals=l[5])
        except:
            bases = "DEL"

        if bases == "DEL": ##DELETION DONT ADD TO STRING
            pass
        elif len(set(bases)) == 1 and  bases[0] == ".": #ALL SAME AS REF
            reference_consensus[ref_name] += ref

        elif len(set(bases)) == 1: ##ALL SAME NUCLEO
            reference_consensus[ref_name] += bases[0]

        elif Counter(bases).most_common()[0][0] == ".":
            reference_consensus[ref_name] += ref

        elif Counter(bases).most_common()[0][0] != "." and len(Counter(bases).most_common()[0][0]) == 1 :
            reference_consensus[ref_name] += Counter(bases).most_common()[0][0]

        ##still need one more if else that handles insertions #and maybe deletions


        else:
            sys.stderr.write(set(bases).difference({"."}))
            sys.stderr.write(Counter(bases).most_common()[0][0])
            sys.stderr.write("{0}\t{1}\t{2}".format(ref_name, pos, bases))
            sys.stderr.write("look if further case should be included")
            sys.exit(1)


for k,v in reference_consensus.iteritems():
    sys.stdout.write(">{0}\n{1}\n".format(k,v))