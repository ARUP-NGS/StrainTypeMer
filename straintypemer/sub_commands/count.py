import gzip
from straintypemer.sub_commands.compare import *
from straintypemer import _ROOT
try:
    import cPickle as pickle
except:
    import pickle


def compare_ard(so, kmer_size=31):
    from Bio import SeqIO
    import jellyfish
    _p = "/home/ksimmon/reference/ard/"
    _p = "/Users/ksimmon/Box Sync/ARUP/strainTypeMer_resources/ard/"
    sys.stderr.write("Retrieving antibiotic resistance genes\n")

    descriptions = {}
    for i in open(_p + "categories.txt"):
        v = i.strip().split("\t")
        name = ".".join(v[0].split(".")[:-1])
        descriptions.update({name: v})

    aro_tags = {}
    for i in open(_p + "AROtags.txt"):
        v = i.strip().split("\t")
        # print v
        aro_tags.update({v[2]: v[1]})
    count = 0

    num_of_sequences = len([i.name for i in SeqIO.parse(_p + "ARmeta-genes.fa", "fasta")])
    for s in SeqIO.parse(_p + "ARmeta-genes.fa", "fasta"):
        count += 1

        sys.stderr.write("\rAnalyzed {0} of {1} antibiotic resistant genes".format(count, num_of_sequences))
        if count != num_of_sequences:
            sys.stderr.flush()
        else:
            sys.stderr.write("\n")

        id = s.description.split(" ")[0]
        species = s.description[s.description.rfind("[") + 1:s.description.rfind("]")]
        aro_tag = [i.split(" ")[0] for i in s.description.split(". ") if "ARO:" in i and "ARO:1000001" not in i]
        # print id, species, descriptions[id][1], ",".join([aro_tags[tag] for tag in aro_tag])
        for j in range(0, len(s.seq) - kmer_size + 1):
            kmer = s.seq[j: j + kmer_size]
            mer = jellyfish.MerDNA(str(kmer))
            mer.canonicalize()

            if id in so.ard:
                so.ard[id][0].append(so.qf[mer])
            else:
                so.ard.update({id: ([so.qf[mer]], species, descriptions[id][1],
                                    [aro_tags[tag] for tag in aro_tag])})
    return

def determine_file_type(this_file, gzipped=False):
    """
    returns 'fa' 'fq' or raises error
    """
    is_type = None
    # get the first line
    if gzipped:
        with gzip.open(this_file, 'r') as f:
            first_line = f.readline()
    else:
        with open(this_file, 'r') as f:
            first_line = f.readline()

    # infer format
    if first_line[0] == "@":
        is_type = "fq"
    elif first_line[0] == ">":
        is_type = "fa"

    # raise error if does not look like valid seq format
    if is_type not in ("fq", "fa"):
        raise TypeError("The files do not appear to be valid 'fasta' of 'fastq' format")
    return is_type


def count_kmers(files, label, gzipped, cpus=1, qual_filter=0, hash_size="500M", no_kmer_filtering=False):
    """
    makes system call to jellyfish to count the kmers in the fastq set

    :param files_to_compare: tuple with (fastq1, fastq2, label, filtering_flag)
    :param gzipped: True if files are gzipped
    :param cpus: number of cpus to pass to jellyfish
    :param qual_filter: flag to pass to jf to filter bases
    :param hash_size: the initial hash size [default "500M"]
    :return: tuple with (label, jf_file_path, filtering_flag)
    """
    _results = []  # will hold a tuple with the 'label', 'jellyfish count file'

    if label is None:
        label = files[0].split(".")[0]

    sys.stderr.write("counting kmers in files {0}: label {1}\n".format(", ".join(files), label))
    # create tmp file
    jf_file = "/tmp/tmp_{0}_{1}.jf".format(label, ''.join(random.choice(string.ascii_uppercase)
                                                          for i in range(8)))

    file_type = determine_file_type(files[0], gzipped=gzipped)
    if no_kmer_filtering:
        count_cutoff = 0
    else:
        count_cutoff = 3

    qual = str(chr(qual_filter + 33))  # convert the qual score into char

    if gzipped:
        if file_type == "fa":  # do not include quality filter
            subprocess.check_call(["gzip -dc {0} | jellyfish count -L {5} -m 31 -s {4} -t {3} -C -o "
                                   "{2} /dev/fd/0".format(" ".join(files), '"' + qual + '"', jf_file, cpus, hash_size,
                                                          count_cutoff)], shell=True)
        elif file_type == "fq":
            subprocess.check_call(["gzip -dc {0} | jellyfish count -Q {2} -L {6} -m 31 -s {5} -t {4} -C -o {3} "
                                   "/dev/fd/0".format(" ".join(files), "", '"' + qual + '"', jf_file, cpus, hash_size,
                                                      count_cutoff)], shell=True)
    else:
        if file_type == "fa":
            subprocess.check_call(["jellyfish", "count", "-L", count_cutoff, "-m", "31", "-s", hash_size,
                                   "-t", str(cpus), "-C", "-o", jf_file] + files)

        if file_type == "fq":
            subprocess.check_call(["jellyfish", "count", "-Q", qual, "-L", count_cutoff, "-m", "31", "-s", hash_size,
                                   "-t", str(cpus), "-C", "-o", jf_file] + files)
    return jf_file


def count(fq_files, gzipped=False, cpus=1, coverage_cutoff=0.15, qual_filter=0, no_kmer_filtering=False,
          label=None, out=sys.stdout):
    """
        This command counts the kmers in a strain and return a comparable list of kmers and other quality
        attributes of the strain.
        :param fq_files: the file list from the commandline (or a string list, see parse_files())
        :param gzipped: True if files are gzipped
        :param cpus: number of cpus to use
        :param coverage_cutoff: percent of genome coverage to set kmer filters [DEFAULT .20]
        :param qual_filter: phred score to filter bases in fq files
        :return: None
    """
    # count kmers in fastq files #these are the raw counts prior to filtering
    # create strain object

    file_path = count_kmers(fq_files, label, gzipped, cpus=cpus, qual_filter=qual_filter)

    so = StrainObject(label, file_path)
    so.do_not_filter = no_kmer_filtering
    so.set_cutoff(coverage_cutoff)

    # Load mlst profiles
    mlst_path = os.path.join(_ROOT, "data/mlst_resources/mlst_profiles.pkl")
    mlst_profiles = None
    if os.path.isfile(mlst_path):
        mlst_profiles = pickle.load(open(mlst_path))

    

    so.filter()
    so.get_stats()

    # antibiotic_resistance_genes
    # if include_ard_comparsion:
    compare_ard(so)

    # TODO write information out to file.

    # so.dump_to_mongo("/home/ksimmon/test_dump.json")
    so.dump_to_mongo(out)

    so.clean_tmp_files()
