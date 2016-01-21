import jellyfish
import subprocess
import sys
import numpy as np
import os
import os.path
from collections import OrderedDict
import itertools

class jf_object:
    """
    Contains information about the jellyfish db for a single strain
    """
    #kmer_cutoff = 0
    jellyfish_path = "jellyfish"

    def __init__(self, name, path):
        """
        initialize a jelly

        :param name:
        :param path:
        :return:
        """
        self.jellyfish_path = "jellyfish"


        self.name = name
        self.path = path
        self.qf = jellyfish.QueryMerFile(path)

        #self.ardb_info_parsed = self.__parser_ardb_info()
        self.get_stats()
        self.histo = self.get_histo()
        self.coverage = self.get_estimate_coverage()
        self.__check_resources()
        self.shared_count = None
        self.kmer_cutoff = None
        self.kmer_reference_count = 0
        self.reverse_kmer_reference_count = 0
        self.kmer_count = 0

        #self.percentage_of_kmers()

    def __check_resources(self):


        if os.path.exists(self.jellyfish_path) is False and os.access(self.jellyfish_path, os.X_OK):
            sys.stderr.write("Jellyfish path not set\n")
            sys.stderr.write(self.jellyfish_path + "\n")
            sys.exit("2")

    def get_stats(self):
        op = subprocess.Popen([self.jellyfish_path, "stats", self.path],
                              stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        out, err = op.communicate()
        if err != "":
            sys.stderr.write("Failed to retrieve db stats\n JELLYFISH ERROR: {0}\n".format(err))
            sys.exit(1)



        out = out.split("\n")
        self.unique_kmers = int(out[0].split(" ")[-1])
        self.distinct_kmers = int(out[1].split(" ")[-1])
        self.total_kmers = int(out[2].split(" ")[-1])
        self.max_count = int(out[3].split(" ")[-1])

    def get_histo(self):

        #print self.max_count
        _histo = OrderedDict()
        #_histo = { i+1 : 0  for i in range(300)}
        op = subprocess.Popen([self.jellyfish_path, "histo", self.path], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        out, err = op.communicate()

        if err != "":
            sys.stderr.write("Failed to retrieve db histo\n JELLYFISH ERROR: {0}\n".format(err))
            sys.exit(1)

        for v in out.strip().split("\n"):
            freq, count = [int(i) for i in v.split(" ")]
            if freq in _histo:
                _histo[freq] += count
            else:
                _histo.update({freq : count })
        return _histo



    def get_estimate_coverage(self):
        """Estimates the coverage for kmers count > 3 times"""
        try:
            if np.sum([float(i) * self.histo[i] for i in self.histo if i >= 3]) == 0 or \
                            float(np.sum(self.histo.values()[3:])) == 0:
                _cov = 1
            else:
                _cov = np.sum([float(i) * self.histo[i] for i in self.histo if i >= 3]) / float(np.sum(self.histo.values()[3:]))
        except:
            _cov = 1
        return _cov

    def estimate_genome_size(self, coverage_cutoff):
        """
        The estimate genome size based on kmer content
        :param coverage_cutoff: the count at which kmers are excluded
        :return: the number of distinct kmers in the strain above cutoff
        """
        return  (np.sum(self.histo.values()[int(coverage_cutoff):]))


    def set_cutoff(self, cutoff):
        self.kmer_cutoff =  cutoff


    def get_shared_count(self):
        return self.shared_count

    def __str__(self):
        return "Name:\t{0}\n".format(self.name) + \
               "Path:\t{0}\n".format(self.path) + \
               "Unique_kmers\t{0}\n".format(self.unique_kmers) + \
               "Distinct_kmers\t{0}\n".format(self.distinct_kmers) + \
               "Total_kmers\t{0}\n".format(self.total_kmers) + \
               "Max_Count\t{0}\n".format(self.max_count)

    def get_kmer_count(self, jellyfish_obj, kmer_reference, reverse_kmer_reference, break_point):
        counter = 0
        _arr = []
        k_arr = None
        r_arr = None

        self.kmer_reference_count = 0
        self.reverse_kmer_reference_count = 0
        self.kmer_count = 0

        if kmer_reference is not None:
            k_arr = []
            kmer_reference = jellyfish.QueryMerFile(kmer_reference)


        if reverse_kmer_reference is not None:
            r_arr = []
            reverse_kmer_reference = jellyfish.QueryMerFile(reverse_kmer_reference)


        mer_file = jellyfish.ReadMerFile(jellyfish_obj)
        for mer, count in mer_file:
            mer.canonicalize()
            _arr.append(self.qf[mer])

            if kmer_reference is not None:
                if kmer_reference[mer] > 0: # kmer in ref
                    k_arr.append(self.qf[mer])

            if reverse_kmer_reference is not None:
                if reverse_kmer_reference[mer] == 0: # kmer not in ref
                    r_arr.append(self.qf[mer])

            counter += 1
            if break_point is not None and counter >= break_point:
                break

        return (self.name, _arr, k_arr, r_arr)



    def mlst_profile(self, mlst_list, kmer_size=31):
        results = ""
        if len(mlst_list) == 0:
            return "no profiles loaded"
        for strain, strain_profiles in mlst_list.iteritems():
            matching_sequences = OrderedDict()
            for gene, sequences in strain_profiles['sequences'].iteritems():
                matching_sequences.update({gene : []})
                for id, sequence in sequences:
                    #print id
                    for j in range(0, len(sequence) - kmer_size + 1):
                        kmer = sequence[j:j+31]
                        mer = jellyfish.MerDNA(str(kmer))
                        mer.canonicalize()
                        if self.qf[mer] == 0: # break loop if not exact match
                            break
                    else:
                        matching_sequences[gene].append(id)
            st_keys =  [":".join(t) for t in list(itertools.product(*matching_sequences.values()))]

            for k in st_keys:
                if k in strain_profiles['st']:
                    st = strain_profiles['st'][k]
                else:
                    st = 'NONE'

                results += "{0}\tST: {2}\tprofile: {1} [{3}]\n".format(strain, k, st,
                                                                       ":".join(strain_profiles['sequences'].keys()))
        return results


# for fasta in fastas:
#         for i, s in enumerate(SeqIO.parse(base + fasta, "fasta")):
#             for j in range(0, len(s) - kmer_size + 1):
#                 kmer = s.seq[j : j+31]
#                 mer = jellyfish.MerDNA(str(kmer))
#                 mer.canonicalize()
#                 if qf[mer] == 0: # break loop if not exact match
#                     break
#             else:
#                 gene = "_".join(s.name.split("_")[0:-1])
#                 seq_num = s.name.split("_")[-1]
#                 matching_sequences[gene].append(seq_num)

#    return [":".join(t) for t in list(itertools.product(*matching_sequences.values()))]


    # def compare_to_ardb(self):
    #     _dict = {}
    #     for _file in os.listdir(self.ardb_dir):
    #         gene = _file.split(".")[0]
    #         for i in open(os.path.join(self.ardb_dir, _file), "r"):
    #             mer, count = i.strip().split("\t")
    #             mer = jellyfish.MerDNA(mer)
    #             mer.canonicalize()
    #             if self.qf[mer] != 0:
    #                 if gene in _dict:
    #                     _dict[gene] += 1
    #                 else:
    #                     _dict.update({gene:1})
    #     return _dict

    # def get_ardb_info(self, k):
    #    return self.ardb_info_parsed[k]
