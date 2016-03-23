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
from itertools import groupby
from collections import Counter

try:
    import jellyfish
except ImportError:
    ImportError.message += "Jellyfish is not installed correctly make sure the python bindings are installed\n"
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
        self.do_not_filter = False
        self.histo = self.get_histo()
        self.coverage = self.get_estimate_coverage()
        self.__check_resources()
        self.shared_count = None
        self.kmer_cutoff = None
        self.kmer_reference_count = 0
        self.estimated_coverage = 0
        self.reverse_kmer_reference_count = 0
        self.kmer_count = 0
        self.kmer_set = set([])
        self.kmer_archive = set([])
        self.filtered_jf_file = "/tmp/tmp_filtered_{0}_{1}.jf".format(self.name ,
                                                                  ''.join(random.choice(string.ascii_uppercase)
                                                                          for i in range(8)))
        self.ard = {}
        self.unique_kmers = None
        self.distinct_kmers = None
        self.total_kmers = None
        self.max_count = None
        self.qf = None
        self.rf = None

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
                _histo.update({freq: count})
        return _histo

    def get_estimate_coverage(self):
        """Estimates the coverage for kmers count > 3 times"""
        try:
            if np.sum([float(i) * self.histo[i] for i in self.histo if i >= 3]) == 0 or \
                            float(np.sum(self.histo.values()[3:])) == 0:
                _cov = 1
            else:
                _cov = np.sum([float(i) * self.histo[i] for i in self.histo if i >= 3]) / float(np.sum(
                        self.histo.values()[3:]))
        except:
            _cov = 1
        return _cov

    def estimate_genome_size(self, coverage_cutoff):
        """
        The estimate genome size based on kmer content
        :param coverage_cutoff: the count at which kmers are excluded
        :return: the number of distinct kmers in the strain above cutoff
        """
        return np.sum(self.histo.values()[int(coverage_cutoff):])

    def set_cutoff(self, cutoff):
        """
        Set the cutoff on which to filter the kmer set
        :param cutoff: int with cutoff value
        :return: None
        """
        self.kmer_cutoff = cutoff
        return

    def filter(self):
        """
        filters the kmer set
        created to allow threading of the filtering
        :return: Name of the object, the kmer set set([]), and the path to the filtered file in tmp directory
        """
        self.__filter_jf_file()
        return self.name, self.kmer_set, self.filtered_jf_file

    def __filter_jf_file(self):
        """
        filters the raw kmer count set based on kmer cutoff and set the queryfile and readfile paths
        :return: None
        """
        if self.do_not_filter:
            self.filtered_jf_file = self.path
            self.qf = jellyfish.QueryMerFile(self.path)
            self.rf = jellyfish.ReadMerFile(self.path)
            self.__create_set()
        else:
            dummy_jf_file = pkg_resources.resource_filename('straintypemer', 'data/dummy_A.jf')
            subprocess.check_call(["jellyfish", "merge", "-L", str(int(self.kmer_cutoff) + 1), "-o",
                                   self.filtered_jf_file, self.path, dummy_jf_file],)

            self.qf = jellyfish.QueryMerFile(self.filtered_jf_file)
            self.rf = jellyfish.ReadMerFile(self.filtered_jf_file)
            self.__create_set()
        return

    def set_jf_file(self, path):
        """
        set the path to the query and readmer files

        :param path:
        :return: None
        """
        self.qf = jellyfish.QueryMerFile(path)
        self.rf = jellyfish.ReadMerFile(path)
        return None

    def __str__(self):
        return "Name:\t{0}\n".format(self.name) + \
               "Path:\t{0}\n".format(self.path) + \
               "Unique_kmers\t{0}\n".format(self.unique_kmers) + \
               "Distinct_kmers\t{0}\n".format(self.distinct_kmers) + \
               "Total_kmers\t{0}\n".format(self.total_kmers) + \
               "Max_Count\t{0}\n".format(self.max_count)

    # def kmer_file(self):
    #     temp_file = "/tmp/tmp_kmer_{0}".format(''.join(random.choice(string.ascii_uppercase) for i in range(8)))
    #     jf_temp_file = temp_file + ".jf"
    #     try:
    #         cmd = "jellyfish dump -L 11 A1A*all*.jf | jellyfish count -m 31 -s 2G -t 6 /dev/fd/0 -o " + jf_temp_file
    #         subprocess.Popen(cmd, shell=True)
    #     except:
    #         sys.stderr.write("Error in running jellyfish count\n")
    #         raise RuntimeError

    def __create_set(self):
        for mer, count in self.rf:
            self.kmer_set.add(str(mer))
            #the repo strain
            self.kmer_archive.add(str(mer))

    def compare_to(self, strain, reference_set=None, inverse=False):
        if reference_set is None:
            strain_1 = self.kmer_set
            strain_2 = strain.kmer_set
        else:
            if inverse:
                strain_1 = self.kmer_set.difference(reference_set)
                strain_2 = strain.kmer_set.difference(reference_set)
            else:
                strain_1 = self.kmer_set.intersection(reference_set)
                strain_2 = strain.kmer_set.intersection(reference_set)


        intersection = float(len(strain_1.intersection(strain_2)))
        denom = ((len(strain_1) - intersection) + (len(strain_2) - intersection)) + intersection

        total = intersection / denom * 100.0
        smallest_count = float(len(strain_1))
        if len(strain_2) < smallest_count:
            smallest_count = len(strain_2)

        rescue = (len(strain_1.intersection(strain_2)) * 2.0) / (smallest_count * 2.0) * 100.0
        return self.name, strain.name, total, rescue, denom, smallest_count

    def compare_to_and_filter(self, strain, complexity_cutoff=12, coverage_cutoff=3,
                              reference_set=None, inverse=False, filtering_cutoff=90):


        if reference_set is None:
            strain_1 = self.kmer_set
            strain_2 = strain.kmer_set
        else:
            if inverse:
                strain_1 = self.kmer_set.difference(reference_set)
                strain_2 = strain.kmer_set.difference(reference_set)
            else:
                strain_1 = self.kmer_set.intersection(reference_set)
                strain_2 = strain.kmer_set.intersection(reference_set)


        intersection = float(len(strain_1.intersection(strain_2)))
        denom = ((len(strain_1) - intersection) + (len(strain_2) - intersection)) + intersection

        total = intersection / denom * 100.0
        smallest_count = float(len(strain_1))
        strain_1_smallest = True
        if len(strain_2) < smallest_count:
            smallest_count = len(strain_2)
            strain_1_smallest = False

        rescue_numerator = float(len(strain_1.intersection(strain_2)))
        rescue = rescue_numerator / (smallest_count) * 100.0
        #return self.name, strain.name, total, rescue, denom, smallest_count
        if total < filtering_cutoff:
            return self.name, strain.name, total, rescue, denom, smallest_count


        # get the difference kmers
        differences_1 = strain_1.difference(strain_2)
        differences_2 = strain_2.difference(strain_1)
        differences = strain_1.symmetric_difference(strain_2)
        # differences = differences_1.union(differences_2) # combined

        complexity_count = 0
        within_1 = 0
        within_2 = 0
        within_3 = 0
        counter_not_filtered = 0
        counter_filtered = 0
        coverage_100 = 0
        kept_1 = 0
        kept_2 = 0
        nucleotide_skew = 0
        filtered_1 = 0
        filtered_2 = 0

        for i, kmer in enumerate(differences):
            is_filtered_kmer = False
            mer = jellyfish.MerDNA(kmer)
            mer.canonicalize()
            s1_count = int(self.qf[mer])
            s2_count = int(strain.qf[mer])

            # s1_out = "{0}:n={1} [cutoff={2}]".format(self.name, s1_count, self.kmer_cutoff)
            # s2_out = "{0}:n={1} [cutoff={2}]".format(strain.name, s2_count, strain.kmer_cutoff)

            if s1_count > self.coverage * coverage_cutoff:
                coverage_100 += 1
                is_filtered_kmer = True
            elif s2_count > strain.coverage * coverage_cutoff:
                coverage_100 += 1
                is_filtered_kmer = True

            if s2_count == 0:
                if  s1_count - int(self.kmer_cutoff)  == 1:
                    within_1 += 1
                    within_2 += 1
                    within_3 += 1
                    is_filtered_kmer = True
                elif s1_count - int(self.kmer_cutoff)  == 2:
                    #within_1 += 1
                    within_2 += 1
                    within_3 += 1
                    is_filtered_kmer = True
                elif s1_count - int(self.kmer_cutoff)  == 3:
                    #within_1 += 1
                    #within_2 += 1
                    within_3 += 1
                    is_filtered_kmer = True
            else:
                if  s2_count - int(strain.kmer_cutoff)  == 1:
                    within_1 += 1
                    within_2 += 1
                    within_3 += 1
                    is_filtered_kmer = True
                elif s2_count - int(strain.kmer_cutoff)  == 2:
                    #within_1 += 1
                    within_2 += 1
                    within_3 += 1
                    is_filtered_kmer = True
                elif s2_count - int(strain.kmer_cutoff)  == 3:
                    #within_1 += 1
                    #within_2 += 1
                    within_3 += 1
                    is_filtered_kmer = True


            complexity = [[k, len(list(g))] for k, g in groupby(kmer)]
            complexity = sorted(complexity,key=lambda l:l[1], reverse=True)
            complexity = sum([v for g, v in complexity[:3]])

            complexity_char = sorted(Counter(kmer).values(), reverse=True)
            if complexity_char[0] > (31.0 / 2):
                nucleotide_skew += 1
                is_filtered_kmer = True

            if complexity > complexity_cutoff:
                is_filtered_kmer = True
                complexity_count += 1

            if is_filtered_kmer:
                counter_filtered += 1
                if s1_count == 0:
                    filtered_2 +=1
                else:
                    filtered_1 +=1
            else: # is_filtered_kmer is False:
                counter_not_filtered += 1
                if s1_count == 0:
                    kept_2 += 1
                else:
                    kept_1 += 1


                #sys.stdout.write(">{0}\t{1}\t{2}\tcomplexity:{3}\n{4}\n".format(
                #    counter_not_filtered, s1_out, s2_out, complexity, kmer))

        print
        print "{0:.1f} coverage {1}".format(self.coverage, self.name)
        print "{0:.1f} coverage {1}".format(strain.coverage, strain.name)
        print "{0} kmers not found in {1} but not {2}".format(len(differences_1), self.name, strain.name)
        print "{0} kmers not found in {1} but not {2} [FILTERED]".format(kept_1, self.name, strain.name)
        print "{0} kmers not found in {1} but not {2}".format(len(differences_2), strain.name, self.name)
        print "{0} kmers not found in {1} but not {2} [FILTERED]".format(kept_2, strain.name, self.name)
        print "3x homopolymer runs summing >= {0} (n={1})".format(complexity_cutoff, complexity_count)
        print "half of kmer one base (n={0})".format(nucleotide_skew)
        print "within 1:2:3 of cutoff= {0}: {1}: {2}".format(within_1, within_2, within_3)
        print "kmers filter = {0}\tkmers retained = {1}".format(counter_filtered, counter_not_filtered)
        print "kmers with excessive coverage [{0}X]= {1}".format(coverage_cutoff ,coverage_100)
        print
        denom -= counter_filtered
        if strain_1_smallest:
            smallest_count -= filtered_1
        else:
            smallest_count -= filtered_2
        total = intersection / denom * 100.0
        rescue = rescue_numerator / (smallest_count) * 100.0
        return self.name, strain.name, total, rescue, denom, smallest_count




    def clean_tmp_files(self):
        try:
            os.remove(self.path)
        except OSError:
            pass
        try:
            os.remove(self.filtered_jf_file)
        except OSError:
            pass

    def mlst_profiles(self, mlst_profiles):
        results = []
        if mlst_profiles is None:
            return ["no profiles loaded"]

        matching_sequences = OrderedDict()
        for species, _d in mlst_profiles.iteritems():
            # print species
            matching_sequences.update({species: {}})
            for i, gene in enumerate(_d["GENE_ORDER"]):
                matching_sequences[species].update({gene: []})
                for profile_number, profile in _d["GENES"][gene].iteritems():
                    if len(profile) == len(profile.intersection(self.kmer_set)):
                        # print 'matched'
                        matching_sequences[species][gene].append(profile_number)

            st_keys = [":".join(t) for t in list(itertools.product(*matching_sequences[species].values()))]
            # print st_keys
            for k in st_keys:
                if k in _d["ST"]:
                    st = _d["ST"][k]
                else:
                    st = 'NONE'
                results.append("{0}\tST: {2}\tprofile: {1} [{3}]".format(species, k, st, ":".join(_d["GENE_ORDER"])))
        return results

    def ard_result(self, coverage_cutoff=0.80):
        _out = {}
        for id, result in self.ard.iteritems():
            counts = np.array(result[0]).clip(0,1)
            gene_length = len(counts)
            kmers_covered = float(np.sum(counts))
            per_covered = kmers_covered / gene_length
            tag = ";".join(result[3])
            if per_covered >= coverage_cutoff:
                if tag not in _out:
                    _out.update({ tag :{"percent_covered": per_covered,
                             "gene_length": gene_length,
                             "ref_id": id,
                             "tag": tag,
                             "description": result[2],
                             "species": result[1] }
                    })
                else:
                    if per_covered > _out[tag]["percent_covered"] and kmers_covered > _out[tag]["gene_length"]:
                        _out[tag]["percent_covered"] = per_covered
                        _out[tag]["gene_length"] = gene_length
                        _out[tag]["ref_id"] = id
                        _out[tag]["species"] = result[1]
        return _out
