#!/home/ksimmon/anaconda/bin/python

__author__ = 'keith simmon'

"""
This script will take as input jf databases from strains and perform a comparison

The input of strains should be provided should be provided as a comma separated list.
assume sample name is before underscore
"""
jellyfish_path = "/home/ksimmon/bin/Jellyfish/bin/jellyfish"
threads = 8


import jellyfish
import numpy as np
import argparse
import os.path
import sys
import subprocess
import collections
import jf_object as jfobj
from multiprocessing import Process, Queue

def count_kmers(q, merged_jf_obj, this_jf_object, cutoff):
    """
    This function facilitates the queuing of the jellyfish counting

    Takes a jellyfish object merged from all the strains and calls the gets kmer count for the queried strain

    :param q: queue
    :param merged_jf_obj: merged count of all strains
    :param this_jf_object: strain jf object
    :return: None (puts in queue)    """
    q.put(this_jf_object.get_kmer_count(merged_jf_obj, cutoff))

def main():
    """
    The main method which takes a list jf counts (strains) compares them.
    :return: status
    """
    parser = argparse.ArgumentParser(description="the script take multiple jellyfish counts and compares the strains")
    parser.add_argument("-c", "--cutoff", help="The number of kmers to anaylaze [Default = None]", type=int, default=None)
    parser.add_argument('jf_files', nargs='+', help='jellyfish files for each strain')
    args = parser.parse_args()
    cutoff = args.cutoff
    jf_files = args.jf_files
    ########################################################################################################################

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

    #######################################################################################################################
    #TODO create unique name that to avoid overwriting when running multiple instances
    #TODO clean these files up
    #create merged file
    try:
        subprocess.check_call([jellyfish_path,
                               "merge", "-o", "-L 2", "/tmp/tmp.jf"] + file_paths )
                               #TODO create unique name in tmp and make sure you clean up after yourself
    except:
        sys.stderr.write("Error in running jellyfish merge\n")
        sys.exit(3)

    ########################################################################################################################
    #START THE WORK
    ########################################################################################################################
    count_table = [[] for i in range(len(strain_objs)) ]
    kmer_table = []
    counter = 0

    attach = count_table.append


    merged_jf = jellyfish.ReadMerFile("/tmp/tmp.jf")

    q = Queue()


    jobs = []
    for _obj in strain_objs.itervalues():
        p = Process(target=count_kmers, args=( q, merged_jf, _obj, cutoff, ))
        jobs.append(p)
    for j in jobs:
        j.start()
    for j in jobs:
        _name, _arr = q.get()
        count_table[strain_objs.keys().index(_name)]=np.array(_arr)
        _arr = ""
    q.close()
    merged_jf = "" #kill any memory this is holding

    lowest_coverage = 10000 #set too high to be realistic
    sys.stdout.write(" COVERAGE INFORMATION ".center(80, "-") + "\n")
    sys.stdout.write("calculations only include kmers counted > 3 times\n")


    #determine kmer cutoff for each strain ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #TODO Throw warning if below X
    cutoff_table = [] #this will hold the cutoff
    for k, v in strain_objs.iteritems():
        sys.stdout.write("Strain: {:s}\t Coverage Estimate: {:.1f}\n".format(k, v.coverage))
        cutoff_table.append(int(v.coverage * .25) + 1)
        v.set_cutoff(int(v.coverage * .25) + 1)
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~





    sys.stdout.write(" STRAIN STATS ".center(80, "-") + "\n")
    for k, v in strain_objs.iteritems():
        sys.stdout.write("Strain: {:s}\n".format(k))
        sys.stdout.write("\tInferred genome size: {:.0f}  [filtering kmers counted <= than {:.0f} times]\n".format(
                          v.estimate_genome_size(v.kmer_cutoff), v.kmer_cutoff ))
    sys.stdout.write("".center(80, "-") + "\n\n")


    #FILTER TABLE CREATED
    cutoff = np.array(cutoff_table)

    ##FILTER COUNT TABLE
    count_table_filtered = []
    total_kmers = 0.0
    kept_kmers = 0.0

    for kmer_count in np.array(count_table).T:
        print kmer_count
        total_kmers += 1
        if (np.array(kmer_count) - cutoff).max() >= 0:
            count_table_filtered.append(kmer_count.clip(0,1))
            kept_kmers += 1

    sys.stdout.write(" FILTER INFO ".center(80, "-") + "\n")
    sys.stdout.write("{:.0f}\tkmers counted (n > 1) in all strains\n".format(total_kmers))
    sys.stdout.write("{:.0f}\tkmers included in analysis\n".format(kept_kmers))
    sys.stdout.write("".center(80, "-") + "\n\n")

    count_table_filtered = np.array(count_table_filtered).T #column become row and rows become columns.
    strain_keys = strain_objs.keys()


    #DETERMINE RELATIONSHIPS~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    similarity_dict = collections.OrderedDict()
    for i in range(len(count_table_filtered)):
        similarity_dict.update({strain_keys[i] : collections.OrderedDict()})
        for j in range(i+1, len(count_table_filtered)):
             sum_1 = sum(count_table_filtered[i])
             sum_2 = sum(count_table_filtered[j])
             intersection = float(sum( (count_table_filtered[i] + count_table_filtered[j]) == 2) )
             total_kmers = sum( (count_table_filtered[i] == count_table_filtered[j]).clip(0,1) )
             smallest = sum_1
             if sum_2 < smallest:
                 smallest = sum_2

             #print intersection, total_kmers
             similarity_dict[strain_keys[i]].update({
                                                     strain_keys[j] :
                                                         (intersection / total_kmers * 100,
                                                          intersection / smallest * 100,)
                                                    })

    #DETERMINE RELATIONSHIPS~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    delimeter = ","
    _str = delimeter + delimeter.join(strain_keys) + "\n"
    for i in range(len(strain_keys)):
        _str += strain_keys[i] + delimeter
        for j in range(len(strain_keys)):
            if i == j:
                _str +=  "100" + delimeter
            elif i < j:
                _str += "{:.1f}{:s}".format(similarity_dict[strain_keys[i]][strain_keys[j]][0],delimeter)
            elif i > j:
                _str += "{:.1f}{:s}".format(similarity_dict[strain_keys[j]][strain_keys[i]][1],delimeter)

        _str  = _str[:-1] + "\n"

    sys.stdout.write(_str + "\n")

if __name__ == "__main__":
    main()