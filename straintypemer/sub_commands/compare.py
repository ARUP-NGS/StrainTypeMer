import os
import sys
import subprocess
import random
import string
import pkg_resources
import gzip
from multiprocessing import Process, Queue
from collections import OrderedDict
from straintypemer.sub_commands.plots import *
from straintypemer.sub_commands import StrainObject
try:
    import cPickle as pickle
except ImportError:
    import pickle


def parse_files(fq_files):
    """
    Parse the file string

    fastq/a files can be input in several ways

    Single file:
    fq
    Single file with label:
    fq:label
    paired_files:
    fq1,fq2
    paired_files with labels:
    fq1,fq2:label
    Additionally [NF] can be added to flag that the sample should not be filtered
    :param fq_files: String that contains file paths and labels
    :return: (tuple with (file_1,file_2,label[NF]))
    """
    _files = []  # this will hold a tuple with (path_fastq_1, path_fastq_2, label)
    for f in fq_files:
        do_not_filter_kmers = False
        if "[NF]" in f:
            do_not_filter_kmers = True
        f = f.replace("[NF]", "")
        label, f1, f2 = None, None, None
        # get the path to first file { must exist }
        f1 = f.split(":")[0].split(",")[0]
        # get the path to the second file { not required }
        if "," in f:
            f2 = f.split(":")[0].split(",")[-1]
        # get the label
        if ":" in f:  # the input includes the label
            label = f.strip().split(":")[-1]
        else:  # label must be retrieved from file name
            if f2 is None:  # only one file
                label = os.path.basename(f1)
            else:  # get
                label = ""
                for s1, s2 in zip(os.path.basename(f1), os.path.basename(f2)):
                    if s1 == s2:
                        label += s1
                    else:
                        break
                if len(label) == 0:
                    raise NameError("The files names do not have consistent names to create labels")
        # check to see if files are valid
        if (os.path.isfile(f1)) is False:
            raise IOError("file is not a valid file path: {0}".format(f1))
        # file 2
        if f2 is not None and (os.path.isfile(f2)) is False:
            raise IOError("file is not a valid file path: {0}".format(f2))
        _files.append((f1, f2, label, do_not_filter_kmers))
    # check to make sure the labels are unique
    if len(_files) != len(set([f[2] for f in _files])):
        raise IOError("file labels are not unique")
    return _files


def determine_file_type(this_file, gzipped=False):
    """
    returns 'fa' 'fq' or raises error
    """
    is_type = None
    # get the first line
    if gzipped:
        with gzip.open(this_file, 'r') as f:
            first_line = f.readline().decode()
    else:
        with open(this_file, 'r') as f:
            first_line = f.readline()  # .decode()

    # what format, peek and find out
    if first_line[0] == "@":
        is_type = "fq"
    elif first_line[0] == ">":
        is_type = "fa"

    # raise error if does not look like valid seq format
    if is_type not in ("fq", "fa"):
        raise TypeError("The files do not appear to be valid 'fasta' of 'fastq' format")
    return is_type


def input_is_jf(files_to_compare):
    """
    Allows the user to use jf files as input.
    Saves time in kmer counting

    currently this is hidden argument

    :param files_to_compare:
    :return:
    """
    _out = []
    for i, files in enumerate(files_to_compare):
        _out.append((files[2], files[0], files[-1],))
    return _out


def count_kmers(files_to_compare, gzipped, cpus=1, qual_filter=0, hash_size="500M", no_kmer_filtering=False):
    """
    makes system call to jellyfish to count the kmers in the fastq set

    :param files_to_compare: tuple with (fastq1, fastq2, label, filtering_flag)
    :param gzipped: True if files are gzipped
    :param cpus: number of cpus to pass to jellyfish
    :param qual_filter: flag to pass to jf to filter bases
    :param hash_size: the initial hash size [default "500M"]
    :param no_kmer_filtering:
    :return: tuple with (label, jf_file_path, filtering_flag)
    """
    _results = []  # will hold a tuple with the 'label', 'jellyfish count file'
    for i, files in enumerate(files_to_compare):
        sys.stderr.write("\rcounting kmers in strain {0} of {1}: {2}{3}".format(i + 1, len(files_to_compare), files[2],
                                                                                " " * 30))
        if i != len(files_to_compare) - 1:
            sys.stderr.flush()
        else:
            sys.stderr.write("\n")
        # create temp file name
        jf_file = "/tmp/tmp_{0}_{1}.jf".format(files[-2], ''.join(random.choice(string.ascii_uppercase)
                                                                  for i in range(8)))
        # if file2 does not exist change to empty string
        file_type = determine_file_type(files[0], gzipped=gzipped)

        if files[1] is None:
            f2 = ""
        else:
            f2 = files[1]

        if files[-1] or no_kmer_filtering:
            count_cutoff = 0
        else:
            count_cutoff = 3

        qual = str(chr(qual_filter + 33))
        if gzipped:
            if file_type == "fa":
                subprocess.check_call(["gzip -dc {0} {1} | jellyfish count -L {6} -m 31 -s {5} -t {4} -C -o "
                                       "{3} /dev/fd/0".format(files[0], f2, '"' + qual + '"', jf_file, cpus, hash_size,
                                                              count_cutoff)], shell=True)
            elif file_type == "fq":
                subprocess.check_call(["gzip -dc {0} {1} | jellyfish count -Q {2} -L {6} -m 31 -s {5} -t {4} -C -o "
                                      "{3} /dev/fd/0".format(files[0], f2, '"' + qual + '"', jf_file, cpus, hash_size,
                                                             count_cutoff)], shell=True)
        else:
            if file_type == "fa":
                subprocess.check_call(["jellyfish", "count", "-L", str(count_cutoff), "-m", "31", "-s", hash_size, "-t",
                                      str(cpus), "-C", "-o", jf_file, files[0], f2],)

            if file_type == "fq":
                subprocess.check_call(
                    ["jellyfish", "count", "-Q", qual, "-L", str(count_cutoff), "-m", "31", "-s", hash_size,
                     "-t", str(cpus), "-C", "-o", jf_file, files[0], f2],)
        _results.append((files[2], jf_file, files[-1]))
    return _results


def compare(fq_files, gzipped=False, cpus=1, coverage_cutoff=0.15, qual_filter=0, output_matrix=True,
            output_histogram=True, output_prefix="", no_kmer_filtering=False, kmer_reference=None,
            inverse_kmer_reference=None, include_ard_comparison=False, pairwise_kmer_filtering=False,
            rapid_mode=True, jf_input=False, clean_files=True):
    """
    The is the entry point for this subcommand:  compares multiple files

    :param fq_files: the file list from the commandline (or a string list, see parse_files())
    :param gzipped: True if files are gzipped
    :param cpus: number of cpus to use
    :param coverage_cutoff: percent of genome coverage to set kmer filters [DEFAULT .20]
    :param qual_filter: phred score to filter bases in fq files
    :param output_matrix: True if matrix pdf is to be output
    :param output_histogram: True if histogram pdf is to be output
    :param output_prefix: The basename to append to histogram and matrix outputs
    :param no_kmer_filtering: True if turning of kmer filtering
                (if '[NF]' in file label the filtering will be turned off)
    :param kmer_reference: NOT IMPLEMENTED YET
    :param inverse_kmer_reference: NOT IMPLEMENTED YET
    :return: None
    """
    # Initial QC of files names
    if len(fq_files) < 2:
        sys.stderr.write("Not enough files to compare\n")
        raise Exception
    files_to_compare = parse_files(fq_files)

    # count kmers in fastq files #these are the raw counts prior to filtering
    if jf_input:
        counts = input_is_jf(files_to_compare)
    else:
        counts = count_kmers(files_to_compare, gzipped, cpus=cpus, qual_filter=qual_filter,
                             no_kmer_filtering=no_kmer_filtering)



    sys.stdout.write("".center(80, "-") + "\n")
    strain_objs = OrderedDict()  # create the strain objects
    # calculated and set the coverage
    for label, file_path, filter_file in counts:
        jf = StrainObject(label, file_path)
        #jf.do_not_filter = filter_file
        strain_objs.update({label: jf})

        if rapid_mode:
            jf.rapid_mode = True

        if no_kmer_filtering or filter_file:
            jf.kmer_cutoff = 0
            jf.do_not_filter = True
            jf.coverage = 1
        else:
            jf.set_cutoff(coverage_cutoff=coverage_cutoff)

        sys.stdout.write("Strain: {0:s}\t Coverage Estimate: {1:.1f}\n".format(label, jf.coverage))
    sys.stdout.write("".center(80, "-") + "\n")
    ##########################################################################

    # Filter the kmers from the strain objs
    strain_objs = filter_coverage(strain_objs, cpus=cpus)

    # Load mlst profiles
    mlst_path = pkg_resources.resource_filename('straintypemer', 'data/mlst_profiles.pkl')
    mlst_profiles = None
    if os.path.isfile(mlst_path):
        mlst_profiles = pickle.load(open(mlst_path, "rb"))


    # antibiotic_resistance_genes
    if include_ard_comparison:
        strain_objs = compare_ard(strain_objs)

    # Print out strain stats
    sys.stdout.write(" STRAIN STATS ".center(80, "-") + "\n")
    for name, file_path, filter_file in counts:
        sys.stdout.write("Strain: {:s}\n".format(name))
        if strain_objs[name].do_not_filter:
            # strain_objs[name].set_cutoff(0)
            sys.stdout.write("\tInferred genome size: {0:,}  [no kmer filtering]\n".format(
                    strain_objs[name].estimate_genome_size(strain_objs[name].kmer_cutoff)))
        else:
            sys.stdout.write("\tInferred genome size: {0:,}  [filtering kmers counted <= {1:.0f} times]\n".format(
                strain_objs[name].estimate_genome_size(strain_objs[name].kmer_cutoff), strain_objs[name].kmer_cutoff))

        for profile in strain_objs[name].mlst_profiles(mlst_profiles):
            sys.stdout.write("\tMLST profile: {0}\n".format(profile))

        for tag, ar_result in strain_objs[name].ard_result(coverage_cutoff=0.8).items():
            sys.stdout.write(
            "\tARD GENE: Gene tag: {0} Covered: {1:.1f}% (size {2}) "
            "\n\t\tcount mean(min|max): {3} ({5}|{6}); {4:.1f}X change from coverage"
            "\n\t\tSpecies: {7} Ref_id: {8}\n\t\tDescription: {9}\n".format(
                    ar_result["tag"],
                    ar_result["percent_covered"] * 100,
                    ar_result["gene_length"],
                    int(ar_result["mean_count"]),
                    ar_result["coverage_difference"],
                    ar_result["min_count"],
                    ar_result["max_count"],
                    ar_result["species"],
                    ar_result["ref_id"],
                    ar_result["description"],))

    sys.stdout.write("".center(80, "-") + "\n")

    # CREATE KMER REFERENCE MATRIX
    if kmer_reference is not None:
        reference_set = load_kmer_reference(kmer_reference)
        matrix_data, cluster_matrix = calculate_matrix(strain_objs, cpus=cpus, reference_set=reference_set,
                                                       inverse=False, pairwise_kmer_filtering=pairwise_kmer_filtering)
        if output_matrix:
            strain_keys = list(strain_objs.keys())
            strain_kmer_counts = {s_key: len(strain_objs[s_key].kmer_set) for i, s_key in enumerate(strain_keys)}


            generage_matrix(strain_keys, strain_keys, cluster_matrix, output_prefix + "_kmer_reference",
                            strain_kmer_counts)
        sys.stdout.write("".center(80, "-") + "\n")

    if inverse_kmer_reference is not None:
        reference_set = load_kmer_reference(inverse_kmer_reference)
        matrix_data, cluster_matrix = calculate_matrix(strain_objs, cpus=cpus, reference_set=reference_set,
                                                       inverse=True, pairwise_kmer_filtering=pairwise_kmer_filtering)
        if output_matrix:
            strain_keys = list(strain_objs.keys())
            strain_kmer_counts = {s_key: len(strain_objs[s_key].kmer_set) for i, s_key in enumerate(strain_keys)}
            generage_matrix(strain_keys, strain_keys, cluster_matrix, output_prefix + "_inverse_kmer_reference",
                            strain_kmer_counts)
        sys.stdout.write("".center(80, "-") + "\n")



    # CREATE DISTANCE MATRIX
    matrix_data, cluster_matrix = calculate_matrix(strain_objs, cpus=cpus,
                                                   pairwise_kmer_filtering=pairwise_kmer_filtering)
    sys.stdout.write("".center(80, "-") + "\n")
    if output_matrix:
        sys.stderr.write("generating_figures\n")
        strain_keys = list(strain_objs.keys())
        strain_kmer_counts = {s_key: len(strain_objs[s_key].kmer_set) for i, s_key in enumerate(strain_keys)}

        generage_matrix(strain_keys, strain_keys, cluster_matrix, output_prefix, strain_kmer_counts)

    if output_histogram:
        new_objs = OrderedDict()
        for k, v in strain_objs.items():
            if v.do_not_filter is False:
                new_objs.update({k:v})

        produce_histograms(new_objs, output_prefix)

    if clean_files:
        for strain in list(strain_objs.values()):
            strain.clean_tmp_files()
        sys.stderr.write("completed analysis\n")


def calculate_matrix(strain_objs, cpus=2, reference_set=None, inverse=False, pairwise_kmer_filtering=False):
    """
    Calculates the matrix by pairwise comparison of strains.

    :param strain_objs: Dictionary of strain objects
    :param cpus: number of processors to use
    :return: matrix cluster with rescue, matrix cluster without rescue
    """
    q = Queue()
    # create the job queue
    jobs = []
    num_of_strains_counted = 0
    current_processes = []
    strain_keys = list(strain_objs.keys())
    comparisons_to_make = 0
    similarity_dict = OrderedDict()
    for i in range(len(strain_keys)):
        similarity_dict.update({strain_keys[i]: {}})
        for j in range(i + 1, len(strain_keys)):
            jobs.append((strain_keys[i], strain_keys[j],))
            comparisons_to_make += 1

    # these set comparisons should thread OK
    if comparisons_to_make < cpus:
        cpus = comparisons_to_make

    for cpu in range(cpus):
        strain_1, strain_2 = jobs.pop()
        p = Process(target=compare_strains, args=(q, strain_objs[strain_1], strain_objs[strain_2], reference_set,
                                                  inverse, pairwise_kmer_filtering), name="{0}:{1}".format(strain_1, strain_2))
        current_processes.append(p)

    # start the jobs for the correct number of cpus
    while len(current_processes) != 0:
        current_processes.pop().start()

    # keep the queue moving
    while len(jobs) != 0:
        strain_1_name, strain_2_name, total, rescue, total_kmers, smallest = q.get()  # PAUSES UNTIL A JOBS RETURN
        similarity_dict[strain_1_name].update({strain_2_name: (total, rescue, total_kmers, smallest)})
        num_of_strains_counted += 1

        sys.stderr.write("\r{0}\tof\t{1}\tcomparisons made {2}:{3}".format(num_of_strains_counted, comparisons_to_make,
                                                                           strain_1_name, strain_2_name))
        sys.stderr.flush()

        # start next job
        strain_1, strain_2 = jobs.pop()
        p = Process(target=compare_strains,
                    args=(q, strain_objs[strain_1], strain_objs[strain_2], reference_set,
                          inverse, pairwise_kmer_filtering),
                    name="{0}:{1}".format(strain_1, strain_2))
        p.start()
    # nothing else to start
    # wait until the queue returns 'ALL THE THINGS'
    while num_of_strains_counted != comparisons_to_make:  # finished processing
        strain_1_name, strain_2_name, total, rescue, total_kmers, smallest = q.get()
        similarity_dict[strain_1_name].update({strain_2_name: (total, rescue, total_kmers, smallest)})  # PAUSES
        num_of_strains_counted += 1
        sys.stderr.write("\r{0}\tof\t{1}\tcomparisons made {2}:{3}".format(num_of_strains_counted, comparisons_to_make,
                                                                           strain_1_name, strain_2_name))
        if num_of_strains_counted != comparisons_to_make:
            sys.stderr.flush()
        else:
            sys.stderr.write("\n")

    sys.stdout.write("".center(80, "-") + "\n")
    q.close()
    

    # PRINT SIMILARITY TABlE~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if reference_set is None:
        sys.stdout.write("[FULL]\n")
    else:
        if inverse:
            sys.stdout.write("[INVERSE KMER REFERENCE]\n")
        else:
            sys.stdout.write("[KMER REFERENCE]\n")


    sys.stdout.write("[SIMILARITY TABLE]\n")
    delimeter = ","
    _str = delimeter + delimeter.join(strain_keys) + "\n"
    matrix_data = []
    matrix_data_cluster = []
    for i in range(len(strain_keys)):
        _str += strain_keys[i] + delimeter
        matrix_data.append([])
        matrix_data_cluster.append([])
        for j in range(len(strain_keys)):
            if i == j:
                matrix_data[i].append(100)
                matrix_data_cluster[i].append(100)
                # _str += "1.000" + delimeter
                _str += "{:.2f}{:s}".format(100, delimeter)
            elif i < j:
                matrix_data[i].append(similarity_dict[strain_keys[i]][strain_keys[j]][0])
                matrix_data_cluster[i].append(similarity_dict[strain_keys[i]][strain_keys[j]][0])
                _str += "{:.2f}{:s}".format(similarity_dict[strain_keys[i]][strain_keys[j]][0], delimeter)
            elif i > j:
                matrix_data[i].append(similarity_dict[strain_keys[j]][strain_keys[i]][1])
                matrix_data_cluster[i].append(similarity_dict[strain_keys[j]][strain_keys[i]][0])

                _str += "{:.2f}{:s}".format(similarity_dict[strain_keys[j]][strain_keys[i]][1], delimeter)
        _str = _str[:-1] + "\n"
    sys.stdout.write(_str)
    sys.stdout.write("[SIMILARITY TABLE END]\n")

    # ADD DENOMINATOR TO OUTPUT
    sys.stdout.write("\n")
    sys.stdout.write("[DENOMINATOR TABLE]\n")
    _str = delimeter + delimeter.join(strain_keys) + "\n"

    for i in range(len(strain_keys)):
        _str += strain_keys[i] + delimeter
        for j in range(len(strain_keys)):
            if i == j:
                if strain_objs[strain_keys[i]].rapid_mode:
                    _str += "{0}{1}".format(len(strain_objs[strain_keys[i]].kmer_archive), delimeter)
                else:
                    _str += "{0}{1}".format(len(strain_objs[strain_keys[i]].kmer_set), delimeter)
            elif i < j:
                _str += "{:.0f}{:s}".format(similarity_dict[strain_keys[i]][strain_keys[j]][2], delimeter)
            elif i > j:
                _str += "{:.0f}{:s}".format(similarity_dict[strain_keys[j]][strain_keys[i]][3], delimeter)
        _str = _str[:-1] + "\n"
    sys.stdout.write(_str)
    sys.stdout.write("[DENOMINATOR TABLE END]\n")
    return matrix_data, matrix_data_cluster


def compare_strains(q, strain_1, strain_2, reference_set=None, inverse=False, pairwise_kmer_filtering=True):
    """
    place strain comparsion in a queue

    :param q: the queue
    :param strain_1:
    :param strain_2:
    :return: None
    """
    if pairwise_kmer_filtering:
        # sys.stdout.write("INCLUDING PAIRWISE FILTER\n")
        q.put(strain_1.compare_to_and_filter(strain_2, reference_set=reference_set, inverse=inverse))
    else:
        q.put(strain_1.compare_to(strain_2, reference_set=reference_set, inverse=inverse))
    return


def filter_coverage(strain_objs, cpus=2):
    """
    Thread function that filters out kmers most likely to be errors, setup jf_obj for comparisons

    :param strain_objs: dictionary of strain objs
    :param cpus:
    :return: Strain objs
    """
    # create the job queue
    q = Queue()
    jobs = []
    current_processes = []
    num_of_strains_filtered = 0  # counter

    strain_keys = list(strain_objs.keys())
    for i in range(len(strain_keys)):
        jobs.append(strain_keys[i])  # holds the label for each strain

    # reset the cpus if they are set too high
    if len(jobs) < cpus:
        cpus = len(jobs)

    # find the initial jobs
    for cpu in range(cpus):
        strain = jobs.pop()
        p = Process(target=filter_strains, args=(q, strain_objs[strain]), name=strain + ":filtering")
        current_processes.append(p)

    # start the jobs for the correct number of cpus
    while len(current_processes) != 0:
        current_processes.pop().start()

    # keep the queue moving
    while len(jobs) != 0:
        name, kmer_set, kmer_archive, jf_path = q.get()  # PAUSES UNTIL A JOBS RETURN
        strain_objs[name].kmer_set = kmer_set
        strain_objs[name].kmer_archive = kmer_archive
        strain_objs[name].set_jf_file(jf_path)
        num_of_strains_filtered += 1
        if strain_objs[name].do_not_filter:
            sys.stderr.write("{0} of {1}\t{2} processed [not filtered]\n".format(
                    num_of_strains_filtered, len(strain_objs),  name,))
        else:
            sys.stderr.write("{0} of {1}\t{2} processed [filtered]\n".format(
                    num_of_strains_filtered, len(strain_objs), name))
        # start next job
        strain = jobs.pop()
        p = Process(target=filter_strains, args=(q, strain_objs[strain]), name=strain + ":filtering")
        p.start()
    # nothing else to start
    # wait until the queue returns 'ALL THE THINGS'
    while num_of_strains_filtered != len(strain_objs):  # finished processing
        name, kmer_set, kmer_archive, jf_path = q.get()  # PAUSES UNTIL A JOBS RETURN
        strain_objs[name].kmer_set = kmer_set
        strain_objs[name].kmer_archive = kmer_archive
        strain_objs[name].set_jf_file(jf_path)

        num_of_strains_filtered += 1
        if strain_objs[name].do_not_filter:
            sys.stderr.write("{0} of {1}\t{2} processed [not filtered]\n".format(
                    num_of_strains_filtered, len(strain_objs),  name,))
        else:
            sys.stderr.write("{0} of {1}\t{2} processed [filtered]\n".format(
                    num_of_strains_filtered, len(strain_objs), name))
    q.close()
    return strain_objs


def filter_strains(q, strain):
    """
    Places filtering into queue

    :param q: the queue
    :param strain: strain object to filter
    :return:
    """
    q.put(strain.filter())
    return


def load_kmer_reference(kmer_reference):
    reference_set = set([])
    for i in open(kmer_reference):
        if i[0] != ">":
            reference_set.add(i.strip())
    return reference_set


def compare_ard(strain_objs, kmer_size=31, coverage_cutoff=.40):
    from Bio import SeqIO
    import jellyfish

    categories_file = pkg_resources.resource_filename('straintypemer', 'data/categories.txt')
    arotags_file = pkg_resources.resource_filename('straintypemer', 'data/AROtags.txt')
    armeta_genes_file = pkg_resources.resource_filename('straintypemer', 'data/ARmeta-genes.fa')
    sys.stderr.write("Retrieving antibiotic resistance genes\n")

    descriptions = {}
    for i in open(categories_file):
        v = i.strip().split("\t")
        name = ".".join(v[0].split(".")[:-1])
        descriptions.update({name : v})

    aro_tags = {}
    for i in open(arotags_file):
        v = i.strip().split("\t")
        #print v
        aro_tags.update({v[2]:v[1]})
    count = 0

    num_of_sequences = len([i.name for i in SeqIO.parse(armeta_genes_file, "fasta")])
    for s in SeqIO.parse(armeta_genes_file, "fasta"):
        count += 1
        sys.stderr.write("\rAnalyzed {0} of {1} antibiotic resistant genes".format(count, num_of_sequences))
        if count != num_of_sequences:
            sys.stderr.flush()
        else:
            sys.stderr.write("\n")

        id =  s.description.split(" ")[0]
        species = s.description[s.description.rfind("[") + 1: s.description.rfind("]")]
        aro_tag = [i.split(" ")[0] for i in s.description.split(". ") if "ARO:" in i and "ARO:1000001" not in i]
        #print id, species, descriptions[id][1], ",".join([aro_tags[tag] for tag in aro_tag])
        for j in range(0, len(s.seq) - kmer_size + 1):
            kmer = s.seq[j : j + kmer_size]
            mer = jellyfish.MerDNA(str(kmer))
            mer.canonicalize()
            for label, so in strain_objs.items():
                if id in so.ard:
                    so.ard[id][0].append(so.qf[mer])
                else:
                    so.ard.update({id : ([so.qf[mer]], species, descriptions[id][1],
                                                [aro_tags[tag] for tag in aro_tag] ) })
    return strain_objs
