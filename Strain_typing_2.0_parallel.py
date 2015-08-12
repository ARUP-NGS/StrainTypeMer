#!/home/ksimmon/anaconda/bin/python

__author__ = 'keith simmon'

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
from multiprocessing import Process, Queue
import string
import random


jellyfish_path = "/home/ksimmon/bin/Jellyfish/bin/jellyfish"
threads = 8

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
    parser.add_argument("-c", "--cutoff", help="The number of kmers to analyze [Default = None]", type=int, default=None)
    parser.add_argument( "--no_kmer_filtering", help="Do not filter kmers based on coverage",
                         action="store_true", default=False)
    parser.add_argument("-k", "--kmer_reference",
                            help="instead of merging strains use kmer reference set for comparison",
                            type=str, default=None)
    parser.add_argument('jf_files', nargs='+', help='jellyfish files for each strain')
    args = parser.parse_args()
    no_kmer_filtering = args.no_kmer_filtering
    cutoff = args.cutoff
    jf_files = args.jf_files
    kmer_reference = args.kmer_reference
    ########################################################################################################################
    if kmer_reference is not None and os.path.isfile(kmer_reference) is False:
        sys.stderr.write("kmer reference file does not exist: {0}\n".format(kmer_reference))
        sys.exit(11)



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
    temp_file =  "/tmp/tmp_{0}".format(''.join(random.choice(string.ascii_uppercase) for i in range(8)))
    jf_temp_file = temp_file + ".jf"
    dump_temp_file = temp_file + ".txt"
    if kmer_reference is None:
    #create merged file
        #merge strain files together


        try:
            subprocess.check_call([jellyfish_path,
                                   "merge", '-L', "2", "-o", jf_temp_file] + file_paths )
        except:
            sys.stderr.write("Error in running jellyfish merge\n")
            sys.exit(3)
        #dump kmers out to a file [Technically I should be able to avoid this however a jellyfish bug causes the loss of
        # canonicalization for merged sets ] dumping out to txt file avoids this.
        try:
            subprocess.check_call([jellyfish_path,
                                   "dump", '-c',  "-t", "-o", dump_temp_file, jf_temp_file] )
        except:
            sys.stderr.write("Error in running jellyfish merge\n")
            sys.exit(4)

        os.remove(jf_temp_file)
        merged_jf = dump_temp_file
    else: #supply my own kmer set.
        merged_jf = kmer_reference
    ########################################################################################################################
    #START THE WORK
    ########################################################################################################################
    count_table = [[] for i in range(len(strain_objs)) ]
    kmer_table = []
    counter = 0
    attach = count_table.append

    q = Queue(maxsize=8)
    jobs = []
    for _obj in strain_objs.itervalues():
        p = Process(target=count_kmers, args=( q, merged_jf, _obj, cutoff, ))
        #print q.qsize()
        jobs.append(p)
    for j in jobs:
        j.start()
    for j in jobs:
        _name, _arr = q.get()
        count_table[strain_objs.keys().index(_name)]=np.array(_arr)
    q.close()

    #cleanup temp files
    if os.path.isfile(dump_temp_file):
        os.remove(dump_temp_file)




    #lowest_coverage = 10000 #set too high to be realistic
    sys.stdout.write(" COVERAGE INFORMATION ".center(80, "-") + "\n")
    sys.stdout.write("calculations only include kmers counted > 3 times\n")


    #determine kmer cutoff for each strain ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #TODO Throw warning if below X
    cutoff_table = [] #this will hold the cutoff
    if no_kmer_filtering:
        for k, v in strain_objs.iteritems():
            sys.stdout.write("Strain: {:s}\t Coverage Estimate: {:.1f}\n".format(k, v.coverage))
            cutoff_table.append(0)
            v.set_cutoff(0)

    else:
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
        #print kmer_count
        total_kmers += 1
        if (kmer_count - cutoff).max() >= 0:
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
        sum_1 = np.sum(count_table_filtered[i])  #sum of kmers strain 1
        for j in range(i+1, len(count_table_filtered)):
             sum_2 = np.sum(count_table_filtered[j]) #sum of kmers strain 2
             intersection = sum( (count_table_filtered[i] + count_table_filtered[j]) == 2)
             total_kmers =  sum( (count_table_filtered[i] + count_table_filtered[j]) > 0 )
             smallest = sum_1
             if sum_2 < smallest:
                 smallest = sum_2
             #print sum_1, sum_2, total_kmers, intersection, smallest
             #print intersection, total_kmers
             similarity_dict[strain_keys[i]].update({
                                                     strain_keys[j] :
                                                         (float(intersection) / total_kmers * 100,
                                                          float(intersection) / smallest * 100,
                                                          total_kmers,
                                                          smallest,
                                                          )
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

    sys.stdout.write("\n\n")
    _str = delimeter + delimeter.join(strain_keys) + "\n"
    for i in range(len(strain_keys)):
        _str += strain_keys[i] + delimeter
        for j in range(len(strain_keys)):
            if i == j:
                _str += "{0}{1}".format(np.sum(count_table_filtered[j]),delimeter)
            elif i < j:
                _str += "{:.1f}{:s}".format(similarity_dict[strain_keys[i]][strain_keys[j]][2],delimeter)
            elif i > j:
                _str += "{:.1f}{:s}".format(similarity_dict[strain_keys[j]][strain_keys[i]][3],delimeter)

        _str  = _str[:-1] + "\n"
    sys.stdout.write(_str + "\n")



if __name__ == "__main__":
    main()