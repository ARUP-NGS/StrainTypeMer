#!/home/ksimmon/anaconda/bin/python

__author__ = 'ksimmon'


"""
This script will take as input jf databases from strains and perform a comparison

The input of strains should be provided should be provided as a comma separated list.
assume sample name is before underscore

"""


import jellyfish
import numpy as np
import argparse
import os.path
import sys
import subprocess
import collections
import jf_object as jfobj
from multiprocessing import Process, Queue #TODO

def count_kmers(strain_name, path_to_kmers, queue): #MULTIPROCESSING
    queue.put(strain_objs[strain_name].kmer_count(strain_name, path_to_kmers))


parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", help="List of jf databases seperated by commas (no_spaces)", type=str,)
parser.add_argument("-c", "--cutoff", help="The number of kmers to anaylaze [Default = None]", type=int, default=None)
args = parser.parse_args()
input_string = args.input

#####################################################
#####TESTING#########################################
#_tmp = "E1_S20_clean_q25_31bp.jf,E2_S21_clean_q25_31bp.jf"#,E3_S22_clean_q25_31bp.jf,E4_S23_clean_q25_31bp.jf,EC_S24_clean_q25_31bp.jf"
#_tmp = "E1_S20_60bp.jf,E2_S21_60bp.jf,E3_S22_60bp.jf,E4_S23_60bp.jf,EC_S24_60bp.jf"
#_tmp = "AC8_S4_clean_q25_31bp.jf,A8_S8_clean_q25_31bp.jf"
#base="/home/ksimmon/data/strain_typing/BactieraRun4_22_2014/enterococcus/jellyfish_31_clean/"
#base = "/home/ksimmon/data/strain_typing/Acineto_combined/jellyfish_31_clean/"

#base="/home/ksimmon/data/strain_typing/BactieraRun4_22_2014/enterococcus/jellyfish_60_clean/"
#_tmp = "S1_S12_clean_q25_31bp.jf,S2_S13_clean_q25_31bp.jf,S3_S14_clean_q25_31bp.jf,S4_S15_clean_q25_31bp.jf,S5_S16_clean_q25_31bp.jf,S6_S17_clean_q25_31bp.jf,S7_S18_clean_q25_31bp.jf,SC_S19_clean_q25_31bp.jf"
#base = "/home/ksimmon/data/strain_typing/BactieraRun4_22_2014/staph/jellyfish_31_clean/"
# input_string = _tmp
# for v in _tmp.split(","):
#     input_string += "{0}{1},".format(base, v)
# input_string = input_string[:]
####################################################
####################################################
#base="/home/ksimmon/data/strain_typing/reference_strains/acinetobacter_genome/"
# genomes = [
#     "Acinetobacter_ADP1_uid61597/NC_005966.jf",
#     "Acinetobacter_baumannii_1656_2_uid158677/NC_017162.jf",
#     "Acinetobacter_baumannii_AB307_0294_uid59271/NC_011595.jf",
#     "Acinetobacter_baumannii_ACICU_uid58765/NC_010611.jf",
#     "Acinetobacter_baumannii_ATCC_17978_uid58731/NC_009085.jf",
#     "Acinetobacter_baumannii_AYE_uid61637/NC_010410.jf",
#     "Acinetobacter_baumannii_BJAB07104_uid210971/NC_021726.jf",
#     "Acinetobacter_baumannii_BJAB0715_uid210972/NC_021733.jf",
#     "Acinetobacter_baumannii_BJAB0868_uid210973/NC_021729.jf",
#     "Acinetobacter_baumannii_D1279779_uid190222/NC_020547.jf",
#     "Acinetobacter_baumannii_MDR_TJ_uid162739/NC_017847.jf",
#     "Acinetobacter_baumannii_MDR_ZJ06_uid158685/NC_017171.jf",
#     "Acinetobacter_baumannii_SDF_uid61601/NC_010400.jf",
#     "Acinetobacter_baumannii_TCDC_AB0715_uid158679/NC_017387.jf",
#     "Acinetobacter_baumannii_TYTH_1_uid176498/NC_018706.jf",
#     "Acinetobacter_baumannii_ZW85_1_uid231518/NC_023028.jf",
#     "Acinetobacter_calcoaceticus_PHEA_2_uid83123/NC_016603.jf",
#     "Acinetobacter_oleivorans_DR1_uid50119/NC_014259.jf",
# ]
# for v in genomes:
#     input_string += "{0}{1},".format(base, v)
# input_string = input_string[:-1]
# print input_string

#sys.exit(1)
NUM_OF_STRAINS = len(input_string.split(",")) #calculate the number of samples
file_paths = []

strain_objs = collections.OrderedDict()  ##KEEP RECORD OF FILE PATHS and CREATE THE STRAIN OBJS FOR EACH STRAIN
for f in input_string.strip().split(","):
    file_paths.append(f)
    strain_name = "".join(os.path.basename(f).split("_")[0:2])
    strain_objs.update({strain_name : jfobj.jf_object(strain_name, f)})

    if os.path.exists(f) == False: ##MAKE SURE THE FILE EXISTS
        sys.stderr.write("file does not exist: {0}\n".format(f))
        sys.exit(1)


if len(strain_objs) != NUM_OF_STRAINS: ##MAKE SURE THE NAMES ARE UNIQUE
    sys.stderr.write("strain names are not unique:\n{0}\n".format(", ".join(strain_objs.keys())))
    sys.exit(1)

#######################################################################################################################
#TODO create unique name that to avoid overwriting when running multiple instances
#TODO clean these files up
#create merged file
try:
    subprocess.check_call(["/home/ksimmon/bin/jellyfish-2.2.0/bin/jellyfish", "merge", "-o", "/tmp/tmp.jf"] + file_paths )
except:
    pass


#create_dumped_file
try:
    subprocess.check_call(["/home/ksimmon/bin/jellyfish-2.2.0/bin/jellyfish", "dump", "-tc",  "-L", "2", "-o", "/tmp/dump.fa", "/tmp/tmp.jf"] )
except:
    pass#TODO
########################################################################################################################

count_table = []
kmer_table = []
count = 0
#wf = open("/home/ksimmon/data/E1_E2_unique_kmers_31bp.txt", "w") #TODO REMOVE

with open("/tmp/dump.fa") as f:
    for l in f:
        _arr = []
        attach = _arr.append #is this really faster
        count += 1
        kmer = l.strip().split("\t")[0]

        #if kmer[27:-27] == "CCCGGG": #TODO REMOVE
        mer = jellyfish.MerDNA(kmer)
        mer.canonicalize()

        for v in strain_objs.itervalues():
            attach(int(v.qf[mer])) #query the kmer against each strain THIS can be threaded


        kmer_table.append(kmer)
        count_table.append(_arr)

        #TODO DELETE
        #if  np.array(_arr).min() == 0:
        # if np.array(_arr).max() > 4:
        #     print "{0}\t{1}".format(kmer,  "\t".join([str(i) for i in _arr]))
        #         #wf.write("{0}\t{1}\t{2}\n".format(kmer,_arr[0],_arr[1]))

        #TODO DELETE
        if args.cutoff is not None and count >= args.cutoff:
            break

lowest_coverage = 10000 #too high to be realistic

sys.stdout.write(" COVERAGE INFORMATION ".center(80, "-") + "\n")
sys.stdout.write("calculations only include kmers counted > 3 times\n")


#determine kmer cutoff for each strain
#TODO Throw warning if below X
cutoff = [] #this will hold the cutoff
for k, v in strain_objs.iteritems():
    sys.stdout.write("Strain: {:s}\t Coverage Estimate: {:.1f}\n".format(k, v.coverage))
    if ".jf" not in k:
        cutoff.append(int(v.coverage * .20) + 1)
        v.set_cutoff(int(v.coverage * .20) + 1)
    else:
        cutoff.append(0)
        v.set_cutoff(0)




count_table = np.array(count_table)
sys.stdout.write(" STRAIN STATS ".center(80, "-") + "\n")
for k, v in strain_objs.iteritems():
    sys.stdout.write("Strain: {:s}\n".format(k))
    sys.stdout.write("\tInferred genome size: {:.0f}  [filtering kmers counted <= than {:.0f} times]\n".format(
                      v.estimate_genome_size(v.kmer_cutoff), v.kmer_cutoff ))
sys.stdout.write("".center(80, "-") + "\n\n")


##FILTER TABLE CREATED
cutoff = np.array(cutoff)
##FILTER COUNT TABLE
count_table_filtered = []
total_kmers = 0.0
kept_kmers = 0.0
for kmer_count in count_table:
    total_kmers += 1
    if (np.array(kmer_count) - cutoff).max() >= 0:
        count_table_filtered.append(kmer_count)
        kept_kmers += 1


sys.stdout.write(" FILTER INFO ".center(80, "-") + "\n")
sys.stdout.write("{:.0f}\tkmers counted (n > 1) in all strains\n".format(total_kmers))
#sys.stdout.write("{:.0f} ({:.1f}%)\tkmers remaining after filter in all strains\n".format(kept_kmers, float(kept_kmers)/total_kmers * 100.0))
sys.stdout.write("".center(80, "-") + "\n\n")

count_table_filtered = np.array(count_table_filtered).T #column become row and rows become columns.
strain_keys = strain_objs.keys()

#DETERMINE RELATIONSHIPS
similarity_dict = collections.OrderedDict()
for i in range(len(count_table_filtered) ):
    similarity_dict.update({strain_keys[i] : collections.OrderedDict()})
    for j in range( i + 1, len(count_table_filtered)):
        #print strain_keys[j], strain_keys[i]
        similarity_dict[strain_keys[i]].update({strain_keys[j] :
            float(sum(count_table_filtered[i].astype(bool) == count_table_filtered[j].astype(bool))) /
            float(len(count_table_filtered[0])) * 100 })

delimeter = ","

_str = delimeter + delimeter.join(strain_keys) + "\n"
for i in range(len(strain_keys)):
    _str += strain_keys[i] + delimeter
    for j in range(len(strain_keys)):
        if i == j:
            _str +=  "100" + delimeter
        elif i < j:
            _str += "{:.1f}{:s}".format(similarity_dict[strain_keys[i]][strain_keys[j]],delimeter)
        elif i > j:
            _str += "{:.1f}{:s}".format(similarity_dict[strain_keys[j]][strain_keys[i]],delimeter)

    _str  = _str[:-1] + "\n"

sys.stdout.write(_str + "\n")