import math
import os.path
from multiprocessing import Process, Queue
from straintypemer import _ROOT
from straintypemer.sub_commands.plots import *
from straintypemer.sub_commands import *
from Bio import SeqIO


def mlst_profiles():
    mlst_profile_data = {}
    resource_dir = os.path.join(_ROOT, 'mlst_resources')
    strain_library = [i for i in os.listdir(resource_dir) if i != '__init__.py']
    for strain in strain_library:
        strain_dir = os.path.join(resource_dir, strain)
        resource_files = (os.listdir(strain_dir))
        strain_profile_file = [f for f in resource_files if ".txt" in f][0]
        strain_st = {}
        profile_order = []
        for i, line in enumerate(open(os.path.join(strain_dir, strain_profile_file))):
            if i == 0:
                profile_order = [pro.lower() for pro in line.strip().split("\t")[1: len(resource_files)] ]
            else:
                st = line.strip().split("\t")[0: len(resource_files)]
                strain_st.update({':'.join(st[1:]): st[0]})

        mlst_sequences = OrderedDict()
        for profile in profile_order:
            if profile not in mlst_sequences:
                mlst_sequences.update({profile: []})

            fh = SeqIO.parse(os.path.join(strain_dir, profile + ".tfa"), "fasta")
            for sequence in fh:
                mlst_sequences[profile].append( (sequence.name.split("_")[-1].split('-')[-1], str(sequence.seq), ))
        mlst_profile_data.update({strain : {"st" : strain_st, "sequences" : mlst_sequences}})
    return mlst_profile_data

def compare_strains(jf_files, no_kmer_filtering, cutoff, cpus, coverage_cutoff, kmer_reference=None,
                    inverse_kmer_reference=None, output_histogram=True, output_matrix=True,
                    output_prefix=""):
    """
    This in the primary function which compares the strains

    :param jf_files: [list of file paths (str)]  List of the paths of the jf_files to be compared
    :param no_kmer_filtering: [Boolean] Do not filter out kmers (bad idea)
    :param cutoff: [int] number of kmers to compare (None = all)
    :param cpus: [int] number of processors to utilize
    :param coverage_cutoff: [int] User overide of automated cutoff filter based on estimated sequence coverage
                            (not recommended)
    :param kmer_reference: supplement use kmer reference set for comparison (e.g. plasmid, core genome, pangenome kmers)
    :param inverse_kmer_reference: supplement use kmer reference set for comparison
            (e.g. plasmid, core genome, pangenome kmers)
    :param output_histogram: [boolean] produce a PDF of the coverage histograms
    :param output_matrix: [boolean] produce a PDF of the matrix
    :param output_prefix: [str] Prefix to append to PDF filenames
    :return: None
    """
    ####################################################################################################################
    # CHECK ARGUMENTS AND ATTRIBUTES
    ####################################################################################################################
    if kmer_reference is not None:
        # for ref in kmer_reference:
        if os.path.isfile(kmer_reference) is False:
            sys.stderr.write("kmer reference file does not exist: {0}\n".format(kmer_reference))
            sys.exit(11)

    if inverse_kmer_reference is not None:
        # for ref in reverse_kmer_reference:
        if os.path.isfile(inverse_kmer_reference) is False:
            sys.stderr.write("reverse kmer reference file does not exist: {0}\n".format(inverse_kmer_reference))
            sys.exit(11)

    num_of_strains = len(jf_files)  # calculate the number of samples
    if cpus is None or cpus > num_of_strains:
        cpus = num_of_strains

    file_paths = []

    strain_objs = collections.OrderedDict()  # KEEP RECORD OF FILE PATHS and CREATE THE STRAIN OBJS FOR EACH STRAIN
    for f in jf_files:
        if os.path.exists(f) is False:  # MAKE SURE THE FILE EXISTS
            sys.stderr.write("file does not exist: {0}\n".format(f))
            sys.exit(1)
        file_paths.append(f)
        strain_name = "".join(os.path.basename(f).split(".")[0])

        strain_objs.update({strain_name: jf_object(strain_name, f)})

    if len(strain_objs) != num_of_strains:  # MAKE SURE THE NAMES ARE UNIQUE
        sys.stderr.write("strain names are not unique:\n{0}\n".format(", ".join(strain_objs.keys())))
        sys.exit(2)

    # ##################################################################################################################
    temp_file = "/tmp/tmp_{0}".format(''.join(random.choice(string.ascii_uppercase) for i in range(8)))
    jf_temp_file = temp_file + ".jf"
    # dump_temp_file = temp_file + ".txt"
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #  COMPARE STRAINS DIRECTLY
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    try:
        # subprocess.check_call([jellyfish_path, "stats", '-L', "2",] + [file_paths[0]] )
        subprocess.check_call(["jellyfish", "merge", '-L', "2", "-o", jf_temp_file] + file_paths)
    except:
        sys.stderr.write("Error in running jellyfish merge\n")
        sys.exit(3)

    merged_jf = jf_temp_file
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #  USE USER SUPPLIED KMER LIST  [kmer\tcount]
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    # ##################################################################################################################
    # START THE WORK
    # ##################################################################################################################
    count_table = [[] for i in range(len(strain_objs))]
    count_table_kmer_reference = [[] for i in range(len(strain_objs))]
    count_table_reverse_kmer_reference = [[] for i in range(len(strain_objs))]

    # ~~~~~~~~~~~~~ MULTIPROCESS THIS MOFO ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    q = Queue()  # TODO EXPLORE PIPES
    # create the job queue
    jobs = strain_objs.values()
    num_of_strains_counted = 0
    current_processes = []

    for cpu in range(cpus):
        _obj = jobs.pop()
        p = Process(target=count_kmers, args=(q, merged_jf, _obj, kmer_reference, inverse_kmer_reference, cutoff,),
                    name=_obj.name)
        current_processes.append(p)

    # start the jobs for the correct number of cpus
    while len(current_processes) != 0:
        current_processes.pop().start()

    # keep the queue moving
    while len(jobs) != 0:
        _name, _arr, k_arr, r_arr = q.get()  # PAUSES UNTIL A JOBS RETURN  #AM I WAITING FOR THE FIRST OBJECT
        # ardb_results.update({_name : ardb})
        num_of_strains_counted += 1

        sys.stderr.write("{0}\tof\t{1}\tstrains processed\n".format(num_of_strains_counted, len(strain_objs)))

        count_table[strain_objs.keys().index(_name)] = np.array(_arr)
        count_table_kmer_reference[strain_objs.keys().index(_name)] = np.array(k_arr)
        count_table_reverse_kmer_reference[strain_objs.keys().index(_name)] = np.array(r_arr)

        # start next job
        p = Process(target=count_kmers, args=(q, merged_jf, jobs.pop(), kmer_reference, inverse_kmer_reference,
                                              cutoff,), name=_obj.name)
        p.start()
    # nothing else to start
    # wait until the queue returns 'ALL THE THINGS'
    while num_of_strains_counted != len(strain_objs):  # finished processing
        _name, _arr, k_arr, r_arr = q.get()            # PAUSES UNTIL A JOBS RETURN
        # ardb_results.update({_name : ardb})
        num_of_strains_counted += 1

        sys.stderr.write("{0}\tof\t{1}\tstrains processed\n".format(num_of_strains_counted, len(strain_objs)))
        count_table[strain_objs.keys().index(_name)] = np.array(_arr)
        count_table_kmer_reference[strain_objs.keys().index(_name)] = np.array(k_arr)
        count_table_reverse_kmer_reference[strain_objs.keys().index(_name)] = np.array(r_arr)
    q.close()
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    # cleanup temp files
    if os.path.isfile(merged_jf):
        os.remove(merged_jf)

    sys.stdout.write(" COVERAGE INFORMATION ".center(80, "-") + "\n")
    sys.stdout.write("Calculations only include kmers counted > 3 times\n")

    # determine kmer cutoff for each strain ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    cutoff_table = []  # this will hold the cutoff
    if no_kmer_filtering:
        for k, v in strain_objs.iteritems():
            sys.stdout.write("Strain: {:s}\t Coverage Estimate: {:.1f}\n".format(k, v.coverage))
            v.set_cutoff(0)
    else:
        for k, v in strain_objs.iteritems():
            sys.stdout.write("Strain: {:s}\t Coverage Estimate: {:.1f}\n".format(k, v.coverage))
            cutoff_table.append(math.ceil(v.coverage * coverage_cutoff))
            if math.ceil(v.coverage * coverage_cutoff) <= 3:
                sys.stderr.write(" WARNING ".center(80, "-") + "\n")
                sys.stderr.write("Strain: {0} has low coverage\n".format(k))
                sys.stderr.write("Calculated cutoff is {0}\n".format(int(math.ceil(v.coverage * coverage_cutoff))))
                if (math.ceil(v.coverage * coverage_cutoff)) < 3:
                    sys.stderr.write("Changing kmer cutoff to 3\n")
                    v.set_cutoff(3)
                sys.stderr.write("If estimated genome size is lower than expected consider repeating\n")
                sys.stderr.write("".center(80, "-") + "\n\n")
                v.set_cutoff(int(math.ceil(v.coverage * coverage_cutoff)))
            else:
                if "IT" in v.name:
                    v.set_cutoff(int(math.ceil(v.coverage * (coverage_cutoff * 2) )))
                else:
                    v.set_cutoff(int(math.ceil(v.coverage * coverage_cutoff)))
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    mlst = mlst_profiles() # load mlst profiles if in resource directory

    sys.stdout.write('\n')
    sys.stdout.write(" STRAIN STATS ".center(80, "-") + "\n")
    for k, v in strain_objs.iteritems():
        sys.stdout.write("Strain: {:s}\n".format(k))
        sys.stdout.write("\tInferred genome size: {0:,}  [filtering kmers counted <= {1:.0f} times]\n".format(
            v.estimate_genome_size(v.kmer_cutoff), v.kmer_cutoff))

        # iterate through mlst
        for profile in v.mlst_profile(mlst):
            sys.stdout.write("\tMLST PROFILE: {0}\n".format(profile))






    sys.stdout.write("".center(80, "-") + "\n\n")

    # FILTER TABLE CREATED ie [3,4,5] cutoff for each strian
    cutoff_cov = np.array([v.kmer_cutoff for v in strain_objs.values()])

    # print len(count_table), len(count_table[i])
    # ~~~~~~~~~~~~~~
    # subtract coverage value then clip to binary value (0 absent, 1 present)
    count_table = (np.array(count_table).T - cutoff_cov).clip(0, 1)
    count_table = count_table[~np.all(count_table == 0, axis=1)].T
    # get rid of kmers not present in all strains #memory footprint

    sys.stdout.write(" FILTER INFO ".center(80, "-") + "\n")
    sys.stdout.write(
        "{0:,.0f}\tkmers counted (n > 1) in all strains\n".format(len(count_table[i])))
    sys.stdout.write("".center(80, "-") + "\n\n")

    if kmer_reference is not None:
        sys.stdout.write("[KMER REFERENCE={0}]\n".format(kmer_reference).center(80, "-"))
        count_table_kmer_reference = (np.array(count_table_kmer_reference).T - cutoff_cov).clip(0, 1)
        count_table_kmer_reference = count_table_kmer_reference[~np.all(count_table_kmer_reference == 0, axis=1)].T

        matrix_data_kmer_reference = calculate_matrix(strain_objs, count_table_kmer_reference)
        sys.stdout.write("[END KMER REFERENCE={0}]\n".format(kmer_reference).center(80, "-"))

        if output_matrix:
            strain_keys = strain_objs.keys()
            strain_kmer_counts = {s_key: np.sum(count_table_kmer_reference[i]) for i, s_key in enumerate(strain_keys)}
            generage_matrix(strain_keys, strain_keys, matrix_data_kmer_reference, output_prefix + "_kmer_reference",
                            strain_kmer_counts)

    if inverse_kmer_reference is not None:
        sys.stdout.write("[INVERSE KMER REFERENCE={0}]\n".format(inverse_kmer_reference).center(80, "-"))
        count_table_reverse_kmer_reference = (np.array(count_table_reverse_kmer_reference).T - cutoff_cov).clip(0, 1)
        count_table_reverse_kmer_reference = count_table_reverse_kmer_reference[~np.all(
            count_table_reverse_kmer_reference == 0, axis=1)].T

        matrix_data_reverse_kmer_reference = calculate_matrix(strain_objs, count_table_reverse_kmer_reference)
        sys.stdout.write("[END INVERSE KMER REFERENCE={0}]\n".format(inverse_kmer_reference).center(80, "-"))

        if output_matrix:
            strain_keys = strain_objs.keys()
            strain_kmer_counts = {s_key: np.sum(count_table_reverse_kmer_reference[i])
                                  for i, s_key in enumerate(strain_keys)}
            generage_matrix(strain_keys, strain_keys, matrix_data_reverse_kmer_reference,
                            output_prefix + "_reverse_kmer_reference", strain_kmer_counts)

    # ~~~~~~~~~~~~~~~
    # column become row and rows become columns.
    # strain_keys = strain_objs.keys()

    matrix_data = calculate_matrix(strain_objs, count_table)
    if output_matrix:
        strain_keys = strain_objs.keys()
        strain_kmer_counts = {s_key: np.sum(count_table[i]) for i, s_key in enumerate(strain_keys)}

        generage_matrix(strain_keys, strain_keys, matrix_data, output_prefix, strain_kmer_counts)

    # write out Pdfs

    if output_histogram:
        produce_histograms(strain_objs, output_prefix)
    return


def calculate_matrix(strain_objs, count_table):
    similarity_dict = collections.OrderedDict()
    strain_keys = strain_objs.keys()
    for i in range(len(count_table)):
        similarity_dict.update({strain_keys[i]: collections.OrderedDict()})
        sum_1 = np.sum(count_table[i])  # sum of kmers strain 1
        for j in range(i + 1, len(count_table)):
            sum_2 = np.sum(count_table[j])  # sum of kmers strain 2
            intersection = np.sum((count_table[i] + count_table[j]) == 2)
            total_kmers = np.sum((count_table[i] + count_table[j]) > 0)
            smallest = sum_1
            if sum_2 < smallest:
                smallest = sum_2

            similarity_dict[strain_keys[i]].update({strain_keys[j]:
                                                    (float(intersection) / total_kmers * 100,
                                                     float(intersection) / total_kmers * 100,
                                                     total_kmers, smallest)})

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
                _str += "100" + delimeter
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
                _str += "{0}{1}".format(np.sum(count_table[j]), delimeter)
            elif i < j:
                _str += "{:.0f}{:s}".format(similarity_dict[strain_keys[i]][strain_keys[j]][2], delimeter)
            elif i > j:
                _str += "{:.0f}{:s}".format(similarity_dict[strain_keys[j]][strain_keys[i]][3], delimeter)
        _str = _str[:-1] + "\n"
    sys.stdout.write(_str)
    sys.stdout.write("[DENOMINATOR TABLE END]\n")
    return matrix_data


def count_kmers(q, merged_jf_obj, this_jf_object, kmer_reference, reverse_kmer_reference, cutoff):
    """
    This function facilitates the queuing of the jellyfish counting
    Takes a jellyfish object merged from all the strains and calls the gets kmer count for the queried strain
    :param q: queue
    :param merged_jf_obj: merged count of all strains
    :param this_jf_object: strain jf object
    :param kmer_reference: jf_reference
    :param reverse_kmer_reference: jf_reference
    :param cutoff: number of kmers to analyze
    :return: None
    """
    q.put(this_jf_object.get_kmer_count(merged_jf_obj, kmer_reference, reverse_kmer_reference, cutoff))
    return
