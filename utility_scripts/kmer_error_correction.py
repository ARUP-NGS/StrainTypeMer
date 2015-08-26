__author__ = 'ksimmon'



"""
This script will take as input jf databases from strains and perform a comparison

The input of strains should be provided should be provided as a comma separated list.
assume sample name is before underscore
"""
import jf_object as jfobj

import jellyfish
import numpy as np
from numpy import sum

import argparse
import os.path
import sys
import subprocess
import collections
import string
import random
import itertools

jellyfish_path = "/home/ksimmon/bin/Jellyfish/bin/jellyfish"

def main():
    """
    The main method which takes a list jf counts (strains) compares them.
    :return: status
    """
    parser = argparse.ArgumentParser(description="the script take multiple jellyfish counts and compares the strains")
    parser.add_argument('-k', '--kmer_dump', type= argparse.FileType("r"))
    parser.add_argument('jf_files', nargs='+', help='jellyfish files for each strain')
    args = parser.parse_args()
    #no_kmer_filtering = args.no_kmer_filtering
    #cutoff = args.cutoff
    #cpus = args.cpus
    jf_files = args.jf_files
    #kmer_reference = args.kmer_reference
    ####################################################################################################################
    # CHECK ARGUMENTS AND ATTRIBUTES
    ####################################################################################################################

    NUM_OF_STRAINS = len(jf_files) #calculate the number of samples

    file_paths = []

    strain_objs = collections.OrderedDict()  ##KEEP RECORD OF FILE PATHS and CREATE THE STRAIN OBJS FOR EACH STRAIN
    for f in jf_files:
        if os.path.exists(f) == False: ##MAKE SURE THE FILE EXISTS
            sys.stderr.write("file does not exist: {0}\n".format(f))
            sys.exit(1)
        file_paths.append(f)
        strain_name = "".join(os.path.basename(f).split(".")[0])

        strain_objs.update({strain_name : jfobj.jf_object(strain_name, f)})

    if len(strain_objs) != NUM_OF_STRAINS: ##MAKE SURE THE NAMES ARE UNIQUE
        sys.stderr.write("strain names are not unique:\n{0}\n".format(", ".join(strain_objs.keys())))
        sys.exit(2)


