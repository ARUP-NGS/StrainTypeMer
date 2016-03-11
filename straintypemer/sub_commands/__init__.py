import subprocess
import sys
import numpy as np
import os
import os.path
from collections import OrderedDict
import random
import string
import itertools
import pkg_resources


try:
    import jellyfish
except ImportError:
    sys.stderr.write("Jellyfish is not installed correctly make sure the python bindings are installed\n")
    raise ImportError

class jf_object:
    """
    Contains information about the jellyfish db for a single strain
    """
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
        self.histo = self.get_histo()
        self.coverage = self.get_estimate_coverage()
        self.__check_resources()
        self.shared_count = None
        self.kmer_cutoff = None
        self.kmer_reference_count = 0
        self.reverse_kmer_reference_count = 0
        self.kmer_count = 0
        self.kmer_set = set([])
        self.filtered_jf_file = "/tmp/tmp_filtered_{0}.jf".format(''.join(random.choice(string.ascii_uppercase)
                                                                          for i in range(8)))
        self.qf = None # jellyfish.QueryMerFile(self.filtered_jf_file)
        self.rf = None # jellyfish.QueryMerFile(self.filtered_jf_file)


    def __check_resources(self):
        if os.path.exists(self.jellyfish_path) is False and os.access(self.jellyfish_path, os.X_OK):
            sys.stderr.write("Jellyfish path not set\n")
            sys.stderr.write(self.jellyfish_path + "\n")
            raise ImportError

    def get_stats(self):
        op = subprocess.Popen([self.jellyfish_path, "stats", self.filtered_jf_file],
                              stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        out, err = op.communicate()
        if err != "":
            sys.stderr.write("Failed to retrieve db stats\n JELLYFISH ERROR: {0}\n".format(err))
            raise RuntimeError

        out = out.split("\n")
        self.unique_kmers = int(out[0].split(" ")[-1])
        self.distinct_kmers = int(out[1].split(" ")[-1])
        self.total_kmers = int(out[2].split(" ")[-1])
        self.max_count = int(out[3].split(" ")[-1])

    def get_histo(self):
        _histo = OrderedDict()
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
        self.__filter_jf_file()
        self.get_stats()


    def __filter_jf_file(self):
        dummy_jf_file = pkg_resources.resource_filename('straintypemer', 'data/dummy_A.jf')
        p1 = subprocess.check_call(["jellyfish", "merge", "-L", str(int(self.kmer_cutoff) + 1),
                                    "-o", self.filtered_jf_file, self.path, dummy_jf_file],)
        self.qf = jellyfish.QueryMerFile(self.filtered_jf_file)
        self.rf = jellyfish.ReadMerFile(self.filtered_jf_file)
        self.__create_set()
        return

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


    def kmer_file(self):
        temp_file = "/tmp/tmp_kmer_{0}".format(''.join(random.choice(string.ascii_uppercase) for i in range(8)))
        jf_temp_file = temp_file + ".jf"
        try:
            cmd = "jellyfish dump -L 11 A1A*all*.jf | jellyfish count -m 31 -s 2G -t 6 /dev/fd/0 -o " + jf_temp_file
            subprocess.Popen(cmd, shell=True)
        except:
            sys.stderr.write("Error in running jellyfish count\n")
            raise RuntimeError


    def compare_to(self, strain):
        cmp_jf_file = "/tmp/tmp_cmp_{0}.jf".format(''.join(random.choice(string.ascii_uppercase)
                                                                          for i in range(8)))
        p1 = subprocess.check_call(["jellyfish", "merge", "-o", cmp_jf_file,
                                    self.filtered_jf_file, strain.filtered_jf_file],)

        op = subprocess.Popen([self.jellyfish_path, "stats", cmp_jf_file],
                              stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        out, err = op.communicate()
        if err != "":
            sys.stderr.write("Failed to retrieve db stats\n JELLYFISH ERROR: {0}\n".format(err))
            raise RuntimeError
        out = out.split("\n")
        unique_kmers = int(out[0].split(" ")[-1])
        distinct_kmers = int(out[1].split(" ")[-1])
        total_kmers = int(out[2].split(" ")[-1])
        max_count = int(out[3].split(" ")[-1])
        os.remove(cmp_jf_file)
        return (float(self.distinct_kmers) + strain.distinct_kmers) / (distinct_kmers * 2.0) * 100



    def __create_set(self):
        for mer, count in self.rf:
            self.kmer_set.add(str(mer))

    def compare_to_set(self, strain):

        intersection = float(len(self.kmer_set.intersection(strain.kmer_set)))
        denom = ((len(self.kmer_set) - intersection) + (len(strain.kmer_set) - intersection)) + intersection

        total =  intersection / denom * 100.0
        smallest_count = float(len(self.kmer_set))
        if len(strain.kmer_set) < smallest_count:
            smallest_count = len(strain.kmer_set)

        rescue = (len(self.kmer_set.intersection(strain.kmer_set)) * 2.0) / (smallest_count * 2.0) * 100.0
        return (self.name, strain.name, total, rescue, denom , smallest_count)



    def clean_tmp_files(self):
        os.remove(self.filtered_jf_file)




    def mlst_profiles(self, mlst_profiles):
        results = []
        if mlst_profiles is None:
            return ["no profiles loaded"]

        matching_sequences = OrderedDict()
        for species, _d in mlst_profiles.iteritems():
            for i, gene in enumerate(_d["GENE_ORDER"]):
                matching_sequences.update({gene : []})
                for profile_number, profile in _d["GENES"][gene].iteritems():
                    if len(profile) == len(profile.intersection(self.kmer_set)):
                        matching_sequences[gene].append(profile_number)

            st_keys =  [":".join(t) for t in list(itertools.product(*matching_sequences.values()))]
            for k in st_keys:
                if k in _d["ST"]:
                    st = _d["ST"][k]
                else:
                    st = 'NONE'
                results.append("{0}\tST: {2}\tprofile: {1} [{3}]".format(species, k, st, ":".join(_d["GENE_ORDER"])))
        return results