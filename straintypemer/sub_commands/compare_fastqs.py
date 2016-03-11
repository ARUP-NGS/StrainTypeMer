import argparse
import os
import sys
from multiprocessing import Process, Queue
import path
from exceptions import OSError
from straintypemer.sub_commands.plots import *
import subprocess
import random
import string
import math
from straintypemer.sub_commands import jf_object
import pkg_resources
from straintypemer import _ROOT
import cPickle
from collections import OrderedDict


def parse_files(fq_files):
    """
    Parse the file string
    :param fq_files: String that contains file paths and labels
    :return: (tuple with (file_1,file_2,label) [file_2 can be none
    """
    _files = [] # this will hold a tuple with (path_fastq_1, path_fastq_2, label)
    for f in fq_files:
        label, f1, f2 = None, None, None
        # get the path to first file { must exist }
        f1 = f.split(":")[0].split(",")[0]
        # get the path to the second file { not required }
        if "," in f:
            f2 = f.split(":")[0].split(",")[-1]
        # get the label
        if ":" in f: # the input includes the label
            label = f.strip().split(":")[-1]
        else: # label must be retrieved from file name
            if f2 is None: # only one file
                label = os.path.basename(f1)
            else: # get
                label = ""
                for s1, s2 in zip(os.path.basename(f1), os.path.basename(f2)):
                    if s1 == s2:
                        label += s1
                    else:
                        break
                if len(label) == 0:
                    sys.stderr.write("The files names do not have consistent names to create labels\n")
                    raise StandardError
        # check to see if files are valid
        if (os.path.isfile(f1)) is False:
            sys.stderr.write("file is not a valid file path: {0}\n".format(f1))
            raise IOError
        # file 2
        if f2 is not None and (os.path.isfile(f2)) is False:
            sys.stderr.write("file is not a valid file path: {0}\n".format(f2))
            raise IOError
        _files.append((f1,f2,label))
    # check to make sure the labels are unique
    if len(_files) != len(set([f[2] for f in _files])):
        sys.stderr.write("file labels are not unique\n")
        raise IOError
    return _files


def count_kmers(files_to_compare, gzipped, cpus=1, qual_filter=0):
    _results = [] # will hold a tuple with the 'label', 'jellyfish count file'
    for i, files in enumerate(files_to_compare):
        sys.stderr.write("counting kmers in strain {0}: {1}\n".format(i + 1, files[2]))
        # create temp file name
        temp_file = "/tmp/tmp_{0}.jf".format(''.join(random.choice(string.ascii_uppercase) for i in range(8)))
        # if file2 does not exist change to empty string
        if files[1] is None:
            f2 = ""
        else:
            f2 = files[1]

        qual = '"' + str(chr(qual_filter + 33)) + '"'
        if gzipped:
            p1 = subprocess.check_call(["gzip -dc {0} {1} | jellyfish count -Q {2} -L 3 -m 31 -s 500M -t {4} -C -o "
                                   "{3} /dev/fd/0".format(files[0], f2, qual , temp_file, cpus)], shell=True)
        else:
            p1 = subprocess.check_call(["jellyfish", "count", "-Q", qual, "-L", "3", "-m", "31", "-s", "500M", "-t",
                                        str(cpus), "-C", "-o", temp_file, files[0], f2],)

        _results.append((files[2],temp_file))

    return _results




def compare_fastqs(fq_files=[], gzipped=False, cpus=1, coverage_cutoff=0, qual_filter=0, output_matrix=True,
                   output_histogram=True,
                   output_prefix=""):
    if len(fq_files) < 2:
        sys.stderr.write("Not enough files to compare\n")
        raise Exception
    files_to_compare = parse_files(fq_files)
    counts = count_kmers(files_to_compare, gzipped, cpus=cpus, qual_filter=qual_filter)
    labels = [i[0] for i in counts]

    #print out coverage
    ### THREAD THIS SECTION ##########################################
    strain_objs = {}
    for label, file_path in counts:
        jf = jf_object(label, file_path)
        strain_objs.update({ label : jf })
        sys.stdout.write("Strain: {0:s}\t Coverage Estimate: {1:.1f}\n".format(label, jf.coverage))
        if math.ceil(jf.coverage * coverage_cutoff) <= 3:
            sys.stderr.write(" WARNING ".center(80, "-") + "\n")
            sys.stderr.write("Strain: {0} has low coverage\n".format(label))
            sys.stderr.write("Calculated cutoff is {0}\n".format(int(math.ceil(jf.coverage * coverage_cutoff))))
            if (math.ceil(jf.coverage * coverage_cutoff)) < 3:
                sys.stderr.write("Changing kmer cutoff to 3\n")
                jf.set_cutoff(3)
            sys.stderr.write("If estimated genome size is lower than expected consider repeating\n")
            sys.stderr.write("".center(80, "-") + "\n\n")
            jf.set_cutoff(int(math.ceil(jf.coverage * coverage_cutoff)))
        else:
            jf.set_cutoff(int(math.ceil(jf.coverage * coverage_cutoff)))
        # clean up the tmp file
        os.remove(file_path)
    sys.stdout.write('\n')
    ##########################################################################



    mlst_path = os.path.join(_ROOT, "data/mlst_resources/mlst_profiles.pkl")
    mlst_profiles = None
    if os.path.isfile(mlst_path):
        mlst_profiles = cPickle.load(open(mlst_path))


    # GIVE ME SOME STATS
    sys.stdout.write(" STRAIN STATS ".center(80, "-") + "\n")
    for k, v in strain_objs.iteritems():
        sys.stdout.write("Strain: {:s}\n".format(k))
        sys.stdout.write("\tInferred genome size: {0:,}  [filtering kmers counted <= {1:.0f} times]\n".format(
                                                                v.estimate_genome_size(v.kmer_cutoff), v.kmer_cutoff))
        for profile in  v.mlst_profiles(mlst_profiles):
            sys.stdout.write("\tMLST profile: {0}\n".format(profile))





    #CREATE DISTANCE MATRIX
    matrix_data = calculate_matrix(strain_objs, cpus=cpus)
    if output_matrix:
        sys.stderr.write("generating_figures\n")
        strain_keys = strain_objs.keys()
        strain_kmer_counts = {s_key : len(strain_objs[s_key].kmer_set) for i, s_key in enumerate(strain_keys)}
        generage_matrix(strain_keys, strain_keys, matrix_data, output_prefix, strain_kmer_counts)

    # write out Pdfs

    if output_histogram:
        produce_histograms(strain_objs, output_prefix)


    for strain in strain_objs.itervalues():
        strain.clean_tmp_files()

    sys.stderr.write("completed analysis\n")


def calculate_matrix(strain_objs, cpus=2):
    q = Queue()  # TODO EXPLORE PIPES
    # create the job queue
    jobs = []
    num_of_strains_counted = 0
    current_processes = []

    strain_keys = strain_objs.keys()
    comparisons_to_make = 0
    similarity_dict = OrderedDict()
    for i in range(len(strain_keys)):
        similarity_dict.update({strain_keys[i] : {}})
        for j in range(i + 1, len(strain_keys)):
            jobs.append((strain_keys[i], strain_keys[j],))
            comparisons_to_make += 1

    if jobs < cpus:
        cpus = jobs

    for cpu in range(cpus):
        strain_1, strain_2 = jobs.pop()
        p = Process(target=compare_strains, args=(q, strain_objs[strain_1], strain_objs[strain_2]),
                    name=strain_1 + ":" + strain_2)
        current_processes.append(p)

    # start the jobs for the correct number of cpus
    while len(current_processes) != 0:
        current_processes.pop().start()

    # keep the queue moving
    while len(jobs) != 0:
        strain_1_name, strain_2_name, total, rescue, total_kmers, smallest = q.get()  # PAUSES UNTIL A JOBS RETURN
        similarity_dict[strain_1_name].update({strain_2_name: (total, rescue,total_kmers, smallest)})
        num_of_strains_counted += 1

        sys.stderr.write("{0}\tof\t{1}\tcomparisons made {2}:{3}\n".format(num_of_strains_counted, comparisons_to_make,
                                                                           strain_1_name, strain_2_name))

        # start next job
        strain_1, strain_2 = jobs.pop()
        p = Process(target=compare_strains, args=(q, strain_objs[strain_1], strain_objs[strain_2]),
                     name=strain_1 + ":" + strain_2)
        p.start()
    # nothing else to start
    # wait until the queue returns 'ALL THE THINGS'
    while num_of_strains_counted != comparisons_to_make:  # finished processing
        strain_1_name, strain_2_name, total, rescue, total_kmers, smallest = q.get()
        similarity_dict[strain_1_name].update({strain_2_name: (total, rescue,total_kmers, smallest)})# PAUSES
        num_of_strains_counted += 1
        sys.stderr.write("{0}\tof\t{1}\tcomparisons made {2}:{3}\n".format(num_of_strains_counted, comparisons_to_make,
                                                                           strain_1_name, strain_2_name))

    q.close()





    # PRINT SIMILARITY TABlE~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    sys.stdout.write("[SIMILARITY TABLE]\n")
    delimeter = ","
    _str = delimeter + delimeter.join(strain_keys) + "\n"
    matrix_data = []
    for i in range(len(strain_keys)):
        _str += strain_keys[i] + delimeter

        matrix_data.append([])
        for j in range(len(strain_keys)):
            if i == j:
                matrix_data[i].append(100)
                # _str += "1.000" + delimeter
                _str += "{:.13f}{:s}".format(100, delimeter)
            elif i < j:
                matrix_data[i].append(similarity_dict[strain_keys[i]][strain_keys[j]][0])
                _str += "{:.1f}{:s}".format(similarity_dict[strain_keys[i]][strain_keys[j]][0], delimeter)
            elif i > j:
                matrix_data[i].append(similarity_dict[strain_keys[j]][strain_keys[i]][1])
                _str += "{:.1f}{:s}".format(similarity_dict[strain_keys[j]][strain_keys[i]][1], delimeter)
        _str = _str[:-1] + "\n"
    sys.stdout.write(_str)
    sys.stdout.write("[SIMILARITY TABLE END]\n")

    # ADD DENOMINATOR TO OUTPUT
    sys.stdout.write("\n\n")
    sys.stdout.write("[DENOMINATOR TABLE]\n")
    _str = delimeter + delimeter.join(strain_keys) + "\n"

    for i in range(len(strain_keys)):
        _str += strain_keys[i] + delimeter
        for j in range(len(strain_keys)):
            if i == j:
                _str += "{0}{1}".format(len(strain_objs[strain_keys[i]].kmer_set), delimeter)
            elif i < j:
                _str += "{:.0f}{:s}".format(similarity_dict[strain_keys[i]][strain_keys[j]][2], delimeter)
            elif i > j:
                _str += "{:.0f}{:s}".format(similarity_dict[strain_keys[j]][strain_keys[i]][3], delimeter)
        _str = _str[:-1] + "\n"
    sys.stdout.write(_str)
    sys.stdout.write("[DENOMINATOR TABLE END]\n")
    return matrix_data


def compare_strains(q, strain_1, strain_2 ):
    q.put(strain_1.compare_to_set(strain_2))
    return