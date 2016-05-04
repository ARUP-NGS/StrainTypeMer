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
import math

try:
    import jellyfish
except ImportError:
    ImportError.message += "Jellyfish is not installed correctly make sure the python bindings are installed\n"
    raise ImportError


def rc(sequence):
    _d = {"A": "T", "T": "A", "G": "C", "C": "G"}
    return "".join([_d[s] for s in reversed(sequence)])


class StrainObject:
    """
    Contains information about the jellyfish db for a single strain
    """
    jellyfish_path = "jellyfish"

    def __init__(self, name, path, json_dump=None):
        """
        initialize a StrainObject

        :param name:
        :param path:
        :return:
        """
        if json_dump is None:
            self.jellyfish_path = "jellyfish"
            self.name = name
            self.path = path
            self.do_not_filter = False
            self.histo = self.get_histo()
            self.coverage = self.get_estimate_coverage()
            self.__check_resources()
            self.kmer_cutoff = None
            self.has_suitable_coverage = False
            self.kmer_set = set([])
            self.kmer_archive = set([])
            self.filtered_jf_file = "/tmp/tmp_filtered_{0}_{1}.jf".format(self.name,
                                                                          ''.join(
                                                                              random.choice(string.ascii_uppercase) for
                                                                              i in range(8)))
            self.ard = {}
            self.unique_kmers = None
            self.distinct_kmers = None
            self.total_kmers = None
            self.max_count = None
            self.qf = jellyfish.QueryMerFile(self.path)
            self.qf_filtered = None
            self.rf = None
            self.warnings = []
        else:
            self.name = json_dump["strain_name"]
            self.path = json_dump["path_count_file"]
            self.do_not_filter = json_dump["filtered_kmer_set"]
            self.histo = {int(k): int(v) for k, v in json_dump["kmer_count_histogram"].iteritems()}
            self.coverage = json_dump["coverage"]
            self.kmer_cutoff = json_dump["kmer_cutoff"]
            self.has_suitable_coverage = json_dump["has_suitable_coverage"]
            self.kmer_set = set(json_dump["kmer_archive"])
            self.kmer_archive = set(json_dump["kmer_archive"])
            self.unique_kmers = json_dump["unique_kmers"]
            self.distinct_kmers = json_dump["distinct_kmers"]
            self.max_count = json_dump["max_count"]
            self.warnings = json_dump["warnings"]
            self.ard_result = json_dump["ard_results"]

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

    def set_cutoff(self, coverage_cutoff):
        """
        Set the cutoff on which to filter the kmer set
        This will hold a warning if the strain is below coverage
        :param cutoff: int with cutoff value
        :return: None
        """
        if math.ceil(self.coverage * coverage_cutoff) <= 3:
            self.warnings.append("low coverage")
            self.warnings.append("calculated coverage {0}".format(math.ceil(self.coverage * coverage_cutoff)))
            if (math.ceil(self.coverage * coverage_cutoff)) < 3:
                self.warnings.append("Changed kmer cutoff to 3\n")
                self.kmer_cutoff = 3
            self.warnings.append("If estimated genome size is lower than expected consider repeating")
            self.kmer_cutoff = int(math.ceil(self.coverage * coverage_cutoff))
        else:
            self.kmer_cutoff = int(math.ceil(self.coverage * coverage_cutoff))
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
            self.qf_filtered = jellyfish.QueryMerFile(self.path)
            self.rf = jellyfish.ReadMerFile(self.path)
            self.__create_set()
        else:
            dummy_jf_file = pkg_resources.resource_filename('straintypemer', 'data/dummy_A.jf')
            subprocess.check_call(["jellyfish", "merge", "-L", str(int(self.kmer_cutoff) + 1), "-o",
                                   self.filtered_jf_file, self.path, dummy_jf_file], )

            self.qf_filtered = jellyfish.QueryMerFile(self.filtered_jf_file)
            self.rf = jellyfish.ReadMerFile(self.filtered_jf_file)
            self.__create_set()
        return

    def set_jf_file(self, path):
        """
        set the path to the query and readmer files

        :param path:
        :return: None
        """
        self.qf_filtered = jellyfish.QueryMerFile(path)
        self.rf = jellyfish.ReadMerFile(path)
        return None

    def __str__(self):
        return "Name:\t{0}\n".format(self.name) + \
               "Path:\t{0}\n".format(self.path) + \
               "Unique_kmers\t{0}\n".format(self.unique_kmers) + \
               "Distinct_kmers\t{0}\n".format(self.distinct_kmers) + \
               "Total_kmers\t{0}\n".format(self.total_kmers) + \
               "Max_Count\t{0}\n".format(self.max_count)

    def __create_set(self):
        for mer, count in self.rf:
            self.kmer_set.add(str(mer))
            # the repo strain
            if str(mer)[:2] in {"AG"}:
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
                              reference_set=None, inverse=False, filtering_cutoff=85, verbose=False):

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
        rescue = rescue_numerator / smallest_count * 100.0
        # return self.name, strain.name, total, rescue, denom, smallest_count
        if total < filtering_cutoff:
            return self.name, strain.name, total, rescue, denom, smallest_count

        # get the difference kmers
        differences_1 = strain_1.difference(strain_2)
        differences_2 = strain_2.difference(strain_1)
        differences = strain_1.symmetric_difference(strain_2)
        # differences = differences_1.union(differences_2) # combined

        complexity_count = 0
        within_1_strain_1 = 0
        within_2_strain_1 = 0
        within_3_strain_1 = 0
        within_1_strain_2 = 0
        within_2_strain_2 = 0
        within_3_strain_2 = 0
        counter_not_filtered = 0
        counter_filtered = 0
        coverage_100 = 0
        kept_1 = 0
        kept_2 = 0
        nucleotide_skew = 0
        filtered_1 = 0
        filtered_2 = 0
        below_cutoff_1 = 0
        below_cutoff_2 = 0

        for i, kmer in enumerate(differences):
            is_filtered_kmer = False
            mer = jellyfish.MerDNA(kmer)
            mer.canonicalize()
            s1_count = int(self.qf_filtered[mer])
            s2_count = int(strain.qf_filtered[mer])

            # s1_out = "{0}:n={1} [cutoff={2}]".format(self.name, s1_count, self.kmer_cutoff)
            # s2_out = "{0}:n={1} [cutoff={2}]".format(strain.name, s2_count, strain.kmer_cutoff)

            if s1_count > self.coverage * coverage_cutoff:
                coverage_100 += 1
                is_filtered_kmer = True
            elif s2_count > strain.coverage * coverage_cutoff:
                coverage_100 += 1
                is_filtered_kmer = True

            if s2_count == 0:
                if s1_count - int(self.kmer_cutoff) == 1:
                    within_1_strain_1 += 1
                    within_2_strain_1 += 1
                    within_3_strain_1 += 1
                    is_filtered_kmer = True
                elif s1_count - int(self.kmer_cutoff) == 2:
                    # within_1 += 1
                    within_2_strain_1 += 1
                    within_3_strain_1 += 1
                    is_filtered_kmer = True
                elif s1_count - int(self.kmer_cutoff) == 3:
                    # within_1 += 1
                    # within_2 += 1
                    within_3_strain_1 += 1
                    is_filtered_kmer = True
            else:
                if s2_count - int(strain.kmer_cutoff) == 1:
                    within_1_strain_2 += 1
                    within_2_strain_2 += 1
                    within_3_strain_2 += 1
                    is_filtered_kmer = True
                elif s2_count - int(strain.kmer_cutoff) == 2:
                    # within_1 += 1
                    within_2_strain_2 += 1
                    within_3_strain_2 += 1
                    is_filtered_kmer = True
                elif s2_count - int(strain.kmer_cutoff) == 3:
                    # within_1 += 1
                    # within_2 += 1
                    within_3_strain_2 += 1
                    is_filtered_kmer = True

            complexity = [[k, len(list(g))] for k, g in groupby(kmer)]
            complexity = sorted(complexity, key=lambda l: l[1], reverse=True)
            complexity = sum([v for g, v in complexity[:3]])

            complexity_char = sorted(Counter(kmer).values(), reverse=True)
            if complexity_char[0] > (31.0 / 2):
                nucleotide_skew += 1
                is_filtered_kmer = True

            if complexity > complexity_cutoff:
                is_filtered_kmer = True
                complexity_count += 1

            if is_filtered_kmer is False:  # is_filtered_kmer is False:
                # counter_not_filtered += 1
                if s1_count == 0:
                    if self.qf[mer] == 0:
                        kept_2 += 1
                        if verbose:
                            print "strain: 2\tcount: {0}\t{1}|{2}\t{3}".format(s2_count, kmer, rc(kmer), self.qf[mer])

                    else:
                        below_cutoff_1 += 1
                        is_filtered_kmer = True

                else:
                    if strain.qf[mer] == 0:
                        kept_1 += 1
                        if verbose:
                            print "strain: 1\tcount: {0}\t{1}|{2}\t{3}".format(s1_count, kmer, rc(kmer), strain.qf[mer])

                    else:
                        below_cutoff_2 += 1
                        is_filtered_kmer = True

            # count updating
            if is_filtered_kmer:
                counter_filtered += 1
                if s1_count == 0:
                    filtered_2 += 1
                else:
                    filtered_1 += 1
            else:
                counter_not_filtered += 1

                # sys.stdout.write(">{0}\t{1}\t{2}\tcomplexity:{3}\n{4}\n".format(
                #    counter_not_filtered, s1_out, s2_out, complexity, kmer))
        s = "\n"
        s += "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n"
        s += "{0:.1f} X\tcoverage '{1}' [kmer cutoff = {2}]\n".format(self.coverage, self.name, self.kmer_cutoff)
        s += "{0:.1f} X\tcoverage '{1}' [kmer cutoff = {2}]\n".format(strain.coverage, strain.name, strain.kmer_cutoff)
        s += "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n"
        s += "{0}\tkmers found in '{1}' but not '{2}'\n".format(len(differences_1), self.name, strain.name)
        s += "{0}\tkmers found in '{1}' but not '{2}' [AFTER FILTERING]\n".format(kept_1, self.name, strain.name)
        s += "{0}\tkmers found in '{1}' but not '{2}'\n".format(len(differences_2), strain.name, self.name)
        s += "{0}\tkmers found in '{1}' but not '{2}' [AFTER FILTERING]\n".format(kept_2, strain.name, self.name)
        s += "~~~~~~~~~ FILTER ATTRS ~~~~~~~~~~~\n"
        s += "{1}\tHomopolymer runs summing >= {0} [sum of 3 homopolymer runs]\n".format(complexity_cutoff,
                                                                                         complexity_count)
        s += "{0}\tHalf of the kmer contains a single base\n".format(nucleotide_skew)
        s += "{0} : {1} : {2}\twithin 1:2:3 count of cutoff [strain_1] " \
             "(e.g the kmer is near the histogram tail)\n".format(
            within_1_strain_1, within_2_strain_1, within_3_strain_1)
        s += "{0} : {1} : {2}\tWithin 1:2:3 count of cutoff [strain_2] " \
             "(e.g the kmer is near the histogram tail)\n".format(
            within_1_strain_2, within_2_strain_2, within_3_strain_2)
        s += "{1}\tkmers with excessive coverage [{0}X]\n".format(coverage_cutoff, coverage_100)
        s += "{0}\tkmer found below initial cutoff [strain_1]\n".format(below_cutoff_1)
        s += "{0}\tkmer found below initial cutoff [strain_2]\n".format(below_cutoff_2)
        s += "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n"
        s += "{0}\tkmers filtered\n{1}\tkmers retained\n".format(counter_filtered, counter_not_filtered)
        s += "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n"
        s += "\n"

        # if verbose:
        sys.stdout.write(s)

        denom -= counter_filtered
        if strain_1_smallest:
            smallest_count -= filtered_1
        else:
            smallest_count -= filtered_2
        total = intersection / denom * 100.0
        rescue = rescue_numerator / smallest_count * 100.0
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
            counts = np.array(result[0]).clip(0, 1)
            gene_length = len(counts)
            kmers_covered = float(np.sum(counts))
            per_covered = kmers_covered / gene_length
            tag = ";".join(result[3])
            if per_covered >= coverage_cutoff:
                if tag not in _out:
                    _out.update({tag: {"percent_covered": per_covered,
                                       "gene_length": gene_length,
                                       "ref_id": id,
                                       "tag": tag,
                                       "description": result[2],
                                       "species": result[1]}
                                 })
                else:
                    if per_covered > _out[tag]["percent_covered"] and kmers_covered > _out[tag]["gene_length"]:
                        _out[tag]["percent_covered"] = per_covered
                        _out[tag]["gene_length"] = gene_length
                        _out[tag]["ref_id"] = id
                        _out[tag]["species"] = result[1]
        return _out

    def dump_to_mongo(self, path):
        from bson.json_util import dumps
        from bson import json_util
        contents = OrderedDict()
        contents.update({"strain_name": self.name})
        contents.update({"path_count_file": self.path})
        contents.update({"filtered_kmer_set": self.do_not_filter})
        contents.update({"kmer_count_histogram": self.histo})
        contents.update({"coverage": self.coverage})
        contents.update({"kmer_cutoff": self.kmer_cutoff})
        contents.update({"has_suitable_coverage": self.has_suitable_coverage})
        contents.update({"kmer_count": len(self.kmer_set)})
        contents.update({"kmer_archive_count": len(self.kmer_archive)})
        contents.update({"unique_kmers": self.unique_kmers})
        contents.update({"distinct_kmers": self.distinct_kmers})
        contents.update({"max_count": self.max_count})
        contents.update({"warnings": self.warnings})
        contents.update({"ard_results": self.ard_result(coverage_cutoff=0.90)})
        contents.update({"kmer_archive": list(self.kmer_archive)})
        contents.update({"kmers": list(self.kmer_set)})

        wf = open(path, "w")
        # wf.write(json.dumps(contents, sort_keys=False, indent=4, default=json_util.default))
        wf.write(dumps(contents, sort_keys=False, indent=4, default=json_util.default))
        wf.close()
