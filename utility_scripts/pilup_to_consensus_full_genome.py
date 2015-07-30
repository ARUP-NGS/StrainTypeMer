__author__ = 'ksimmon'
"""
This version of pileup_to_consensus.py take input from a bed file and extracts the pilup
regions for each bedfile and writes fasta files for each bed.


"""
from Bio import SeqIO
import sys
import argparse
from collections import Counter
import numpy as np

def rc(sequence):
    """
    Simple reverse complement method
    :param sequence: input sequence
    :return: reverse complement of the sequence
    """
    _base_key = {"A":"T", "G":"C", "C":"G", "T":"A", "N":"N", "a":"t", "g":"c", "c":"g", "t":"a", "n":"n",}
    return "".join([_base_key[i] for i in reversed(sequence)])
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def get_base(bases, quals, ref):
    """
    :param bases:
    :param quals:
    :param ref:
    :return:
    """
    if len(set(bases)) == 1 and  bases[0] == ".": #ALL SAME AS REF
        return ref, sum(np.array(quals))

    elif len(set(bases)) == 1: ##ALL SAME NUCLEO
        return bases[0], sum(np.array(quals))

    elif Counter(bases).most_common()[0][0] == ".":
        _q = []
        for i, base in enumerate(bases):
            if base == ".":
                _q.append(quals[i])
        return ref, sum(np.array(quals))

    elif Counter(bases).most_common()[0][0] != "." and len(Counter(bases).most_common()[0][0]) == 1 :
        _b = Counter(bases).most_common()[0][0]
        _q = []
        for i, base in enumerate(bases):
            if base == _b:
                _q.append(quals[i])

        return _b, sum(np.array(quals))

    ##still need one more if else that handles insertions #and maybe deletions
    else:
        sys.stderr.write(set(bases).difference({"."}))
        sys.stderr.write(Counter(bases).most_common()[0][0])
        sys.stderr.write("{0}\t{1}\t{2}".format(ref_name, pos, bases))
        sys.stderr.write("look if further case should be included")
        sys.exit(1)





def parse_pileup(coverage=0, bases="", quals="", line="", ref=""):
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
    return get_base(_bases, quals, ref)






def main():

    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", help="mpileup file from samtools", type=str,)# required=True)
    parser.add_argument("-b", "--bed_files", help="bed files separated by a space",
                        nargs="+", type=str)

    args = parser.parse_args()
    pu_file = args.input
    bed_files = args.bed_files

    if pu_file == None and bed_files == None:
        pass
        ##TESTING
        # pu_file = "/home/ksimmon/data/strain_typing/Acineto_combined/mlst3/A3_S3_processed/A3_S3.pileup"
        # bed_files = ["/home/ksimmon/reference/mlst/acinetobacter/oxf_acinetobacter_baumannii_NC_021726_sorted.bed",
        #              "/home/ksimmon/reference/mlst/acinetobacter/pas_acinetobacter_baumannii_NC_021726_sorted.bed"]


    for bed in bed_files:
        #Get all the bed coords from the bed files
        coords = {}
        positions_to_query = set([]) ##this set hold positions from -10 to o
        for line in open(bed, "r"):
            start = int(line.split("\t")[1])
            coords.update({start:[ int(i) if i.isdigit() else i for i in line.strip().split("\t") ] })

        pos_query_key = {}
        for i in coords.values():
            _pos = set([i[1] - j for j in range(0,11)])
            pos_query_key.update({i[1] : _pos})
            positions_to_query.update(_pos)

        reference_consensus = {}
        fh = enumerate(open(pu_file, "r"))
        for i, line in fh:
            l = line.strip().split()
            ref_name, ref_pos, ref, num_reads = l[0], int(l[1]), l[2], int(l[3]) ##pileup
            found_coord = False
            if ref_pos in positions_to_query:
                found_coord = True
                #find the bed coordinate line
                for k,v in pos_query_key.iteritems():
                    if ref_pos in v:
                        current_coord = coords[k]
                ref_coord_name, start_coord, end_coord, gene_coord_name, length, strand = current_coord

                _seq = ""


                ####
                #  GENERATE CONSENSUS
                ##
                while found_coord:
                    try:
                        base, qual = parse_pileup(coverage=num_reads, bases=l[4], quals=l[5], ref=ref)
                    except:
                        base = None


                    #write this mother out
                    if ref_pos < end_coord + 10:
                        if base != None:
                            if (ref_pos < start_coord) or (ref_pos > end_coord):
                                _seq += base.lower()
                            else:
                                _seq += base
                    else:
                        if strand == "+":
                            reference_consensus.update({gene_coord_name :_seq })
                            print ">{0}\n{1}\n".format(gene_coord_name,_seq)
                        else:
                            print ">{0}\n{1}\n".format(gene_coord_name,rc(_seq))
                            reference_consensus.update({gene_coord_name : rc(_seq)})
                        found_coord = False

                    i, l = fh.next()
                    l = l.strip().split("\t")
                    ref_name, ref_pos, ref, num_reads = l[0], int(l[1]), l[2], int(l[3]) ##pileup



if __name__ == "__main__":
    main()