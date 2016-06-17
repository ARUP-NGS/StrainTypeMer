import json
from straintypemer.sub_commands.compare import *


files = [
    "E17_QC_1.json",
    "E17_QC_2.json",
    "USA_300_1.json",
    "staph_qc_1.json",
    "staph_qc_2.json",
    #"USA_300_2.json",
    "aci_qc_1.json",
    "aci_qc_2.json",

]

files = "M10A.json M11A.json M12A.json M13A.json M15A.json M16A.json M17A.json M18A.json M19A.json M1A.json M20A.json " \
        "M2A.json M3A.json M6A.json M7A.json M8A.json M9A.json M10B.json M11B.json M12B.json M13B.json M15B.json " \
        "M16B.json M17B.json M18B.json M19B.json M1B.json M20B.json M2B.json M3B.json M6B.json M7B.json M8B.json " \
        "M9B.json"
files = sorted(files.split(' '))
# files =  [ "_full.".join(f.split(".")) for f in files]

base = "/Users/ksimmon/Box Sync/ARUP/Staph/"
kmer_reference = base + "staph_aureus_80.fa"
strain_objs = collections.OrderedDict()
for f in files:
    s = json.load(open(base + f, "r"))
    print s['strain_name']
    strain_objs.update({s['strain_name'] : StrainObject(name=None, path=None, json_dump=s)})
    print strain_objs[s['strain_name']].histo
matrix_data, cluster_matrix = calculate_matrix(strain_objs, cpus=15)


sys.stderr.write("generating_figures\n")
strain_keys = strain_objs.keys()
strain_kmer_counts = {s_key: len(strain_objs[s_key].kmer_set) for i, s_key in enumerate(strain_keys)}
generage_matrix(strain_keys, strain_keys, cluster_matrix, "staph_cookson", strain_kmer_counts)

produce_histograms(strain_objs, "")


reference_set = load_kmer_reference(kmer_reference)
matrix_data, cluster_matrix = calculate_matrix(strain_objs, cpus=5, reference_set=reference_set,
                                                inverse=False, pairwise_kmer_filtering=False)

strain_keys = strain_objs.keys()
strain_kmer_counts = {s_key: len(strain_objs[s_key].kmer_set) for i, s_key in enumerate(strain_keys)}
generage_matrix(strain_keys, strain_keys, cluster_matrix, "staph_cookson_kmer_reference",
                strain_kmer_counts)
sys.stdout.write("".center(80, "-") + "\n")


matrix_data, cluster_matrix = calculate_matrix(strain_objs, cpus=5, reference_set=reference_set,
                                                inverse=True, pairwise_kmer_filtering=False)

strain_keys = strain_objs.keys()
strain_kmer_counts = {s_key: len(strain_objs[s_key].kmer_set) for i, s_key in enumerate(strain_keys)}
generage_matrix(strain_keys, strain_keys, cluster_matrix, "staph_cookson_inverse_kmer_reference",
                strain_kmer_counts)
sys.stdout.write("".center(80, "-") + "\n")