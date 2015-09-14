__author__ = 'Keith E Simmon'

import os, sys
import os.path
from Bio import SeqIO

def main():

    base_dir = "/home/ksimmon/reference/strian_typing_resources/ARDB/"



    grouped_sequences = {}



    descriptions = {}
    for i in open(base_dir + "class2info.tab"):
        descriptions.update( { i.strip().split("\t")[0] : i.strip().split("\t")[1] } )
        grouped_sequences.update( {i.strip().split("\t")[0] : []} )

    origin_type = {}
    for i in open(base_dir + "origin_type.tab"):
        origin_type.update( { i.strip().split("\t")[0] : i.strip().split("\t")[2] } )

    genes = {}
    for i in open(base_dir + "ar_genes.tab"):
        genes.update( { i.strip().split("\t")[0] : i.strip().split("\t")[1] } )

    sequences = []
    count = 0
    for i in  os.listdir(base_dir + "fasta_sequences/"):

        for s in SeqIO.parse(base_dir + "fasta_sequences/" + i, "fasta"):
            sequences.append( ( i, s,) )


            acc = i.split(":")[1]


            grouped_sequences[origin_type[genes[acc]]].append((i,s,))


    #print count
    #print len(grouped_sequences)
    count = 0
    for k, cl in grouped_sequences.iteritems():
        #print k
        wf = open(base_dir + "grouped_fastas/" + k + ".fa", "w")
        for name, sequence in cl:
            #print name, sequence.name, sequence.seq
            wf.write(">{0}\n{1}\n".format(name[:-3], sequence.seq))
        wf.close()


    #print count

if __name__ == "__main__":
    main()
