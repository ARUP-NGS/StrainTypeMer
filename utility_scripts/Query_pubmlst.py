__author__ = 'keith Simmon'
__date__ = "April 28 2015"
"""
This program takes fasta input of know MLST genes and prints the MLST profile and Stain Type
"""

import sys
from suds.client import Client
import argparse
from StringIO import StringIO


class blast_results():
    """
    a class to hold the blast information from pubMLST
    """
    def __init__(self, locus, id, mismatches, gaps, length, alignment):
        self.locus = locus
        self.id = id
        self.mismatches = mismatches
        self.gaps = gaps
        self.length = length
        self.alignment = alignment

    def is_perfect_match (self):
        if self.gaps == 0 and self.length == self.alignment:
            return True
        else:
            return False

    def locus_to_string(self):
        return "{0} : {1}".format(self.locus, self.id)


##################################################################################################################
def profile_to_string(locus_order, locus_results):
    _s = "["
    for l in locus_order:

        _s += "{0}, ".format(locus_results[l].locus_to_string())
    return _s[:-2] + "]"



def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", help="consensus_fasta_file", type=str,)# required=True)
    parser.add_argument("-d", "--db_to_query", help="the pubmlst db to query", type=str)

    args = parser.parse_args()
    fasta_file = args.input
    database = args.db_to_query

    ## TESTING FILES ##################################################################################################
    if database == None and fasta_file == None:
        pass
        #fasta_file = "/home/ksimmon/data/strain_typing/Acineto_combined/MLST2/A4_S4_consenus.fa"
        #database = "abaumannii"
    ###################################################################################################################

    ## hold fasta sequence as string
    fasta_sequence=""
    for line in open(fasta_file):
        fasta_sequence += line


    #pubmlst priming
    url = "http://pubmlst.org/api/mlst.wsdl"
    client = Client(url, cache=None)

    ##BLAST THE SEQUENCES TO THE DATABASE
    response = client.service.blast(database,  fasta_sequence, 1)
    ###check the results from the blast if not an exact match make a note
    locus_results = {}
    for i in response:
        locus_results.update({i.locus : blast_results(i.locus, i.id, i.mismatches, i.gaps, i.length, i.alignment)})

    ##get the allele list
    locus_list = client.service.getLocusList(database)

    profile = client.factory.create("profile") # create a profile object defined by pubmlst
    for locus in locus_list:
        alleleNumber = client.factory.create("alleleNumber")
        alleleNumber.locus = locus
        alleleNumber.id =  locus_results[locus].id
        profile.alleleNumber.append(alleleNumber)


    sys.stderr = StringIO() # suppress stderr from suds
    try:
        response = client.service.getRelatedSTsByProfile(database, len(locus_list), alleleNumber=profile[0])
        st = response[0]
    except:
        #no ST found that matched the query
        st = None


    sys.stdout.write("ST:\t{0}\t".format(st))
    sys.stdout.write("The MLST Profile:\t{0}\n".format(profile_to_string(locus_list, locus_results)))


if __name__ == '__main__':
    main()