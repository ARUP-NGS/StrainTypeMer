__author__ = 'ksimmon'

import jellyfish
import subprocess
import sys
import numpy as np
import os, os.path
from collections import OrderedDict

class jf_object:
    """
    Contains information about the jellyfish db for a single strain
    """
    kmer_cutoff = 0
    #jellyfish_path = ""# = "/home/ksimmon/bin/jellyfish-2.2.0/bin/jellyfish"
    #ardb_dir = "" # "/home/ksimmon/reference/strian_typing_resources/ARDB/grouped_fastas/jf_files/dump/"
    #ardb_info = ""
    jellyfish_path = "/home/ksimmon/bin/jellyfish-2.2.0/bin/jellyfish"
    #jellyfish_path = jellyfish_path
    ardb_dir = "/home/ksimmon/reference/strian_typing_resources/ARDB/grouped_fastas/jf_files/dump/"
    ardb_info = "/home/ksimmon/reference/strian_typing_resources/ARDB/class2info.tab"


    def __init__(self, name, path):
        """
        initialize a jelly

        :param name:
        :param path:
        :return:
        """
        self.jellyfish_path = "/home/ksimmon/bin/jellyfish-2.2.0/bin/jellyfish"
        self.ardb_dir = "/home/ksimmon/reference/strian_typing_resources/ARDB/grouped_fastas/jf_files/dump/"
        self.ardb_info = "/home/ksimmon/reference/strian_typing_resources/ARDB/class2info.tab"


        self.name = name
        self.path = path
        self.qf = jellyfish.QueryMerFile(path)



        self.ardb_info_parsed = self.__parser_ardb_info()
        self.get_stats()
        self.histo = self.get_histo()
        self.get_estimate_coverage()
        self.__check_resources()
        self.shared_count = None



        #self.percentage_of_kmers()

    def __check_resources(self):
        if os.path.exists(self.jellyfish_path) is False:
            sys.stderr.write("Jellyfish path not set\n")
            sys.stderr.write(self.jellyfish_path + "\n")
            sys.exit("2")
        if os.path.isdir(self.ardb_dir) is False:
            sys.stderr.write("ARDB fasta path not set\n")
            sys.exit("2")
        if os.path.isfile(self.ardb_info)is False:
            sys.stderr.write("ARDB info path not set\n")
            sys.exit("2")

    def __parser_ardb_info(self):
         return {line.split("\t")[0]: line.split("\t")[1].strip() for line in open(self.ardb_info,"r")}

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
            else: #greater than 300
                _histo.update({freq : count })
        return _histo



    def get_estimate_coverage(self):
        """Estimates the coverage for kmers count > 3 times"""
        try:
            #sum of kmers count / the kmers (distinict) count
            self.coverage = np.sum([float(i) * self.histo[i] for i in self.histo]) / float(np.sum(self.histo.values()[3:]))
        except:
            self.coverage = 0


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


    def get_kmer_count(self, jellyfish_obj, break_point):
        counter = 0
        _arr = []

        #add extra info
        ardb = self.compare_to_ardb()


        for i in open(jellyfish_obj, "r"):
            #print i
            mer, count = i.strip().split("\t")
            mer = jellyfish.MerDNA(mer)
            mer.canonicalize()
            #print mer, count, self.qf[mer]
            _arr.append(self.qf[mer])
            counter += 1
            if break_point is not None and counter >= break_point:
                break
        return (self.name, _arr, ardb)

    def compare_to_ardb(self):
        _dict = {}
        for _file in os.listdir(self.ardb_dir):
            gene = _file.split(".")[0]
            for i in open(self.ardb_dir + _file, "r"):
                mer, count = i.strip().split("\t")
                mer = jellyfish.MerDNA(mer)
                mer.canonicalize()
                if  self.qf[mer] != 0:
                    if gene in _dict:
                        _dict[gene] += 1
                    else:
                        _dict.update({gene:1})
        return _dict

    def get_ardb_info(self, k):
        return self.ardb_info_parsed[k]
