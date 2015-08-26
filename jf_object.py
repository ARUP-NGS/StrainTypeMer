__author__ = 'ksimmon'

import jellyfish
import subprocess
import sys
import numpy as np
import multiprocessing

class jf_object:
    """contains information about the jellyfish db for a single strain"""
    kmer_cutoff = 0
    jellyfish_path = "/home/ksimmon/bin/jellyfish-2.2.0/bin/jellyfish"

    def __init__(self, name, path):
        self.name = name
        self.path = path
        self.qf = jellyfish.QueryMerFile(path)

        self.get_stats()
        self.get_histo()
        self.get_estimate_coverage()

        self.shared_count = None
        #self.percentage_of_kmers()


    def get_stats(self):
        op = subprocess.Popen(["/home/ksimmon/bin/jellyfish-2.2.0/bin/jellyfish", "stats", self.path],
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
        op = subprocess.Popen(["/home/ksimmon/bin/jellyfish-2.2.0/bin/jellyfish", "histo", "--threads", "4", "--full", "--high", str(self.max_count ),
                              "--low", "0", self.path], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        out, err = op.communicate()

        if err != "":
            sys.stderr.write("Failed to retrieve db histo\n JELLYFISH ERROR: {0}\n".format(err))
            sys.exit(1)
        _arr = []
        out = out.strip().split("\n")
        self.histo = np.array([ int(i.split(" ")[1])  for i in out])
        self.histo_total = np.array([(i + 1) * self.histo[i] for i in range(len(self.histo))])


    def get_estimate_coverage(self):
        """Estimates the coverage for kmers count > 3 times"""
        try:
            self.coverage = float(np.sum(self.histo_total[3:])) / float(np.sum(self.histo[3:]))
        except:
            self.coverage = 0

    def get_histo_plot(self):
        pass
        #TODO

    def estimate_genome_size(self, coverage_cutoff):  ####GET THE FREQUENCY OF KMERS REPEATED IN ACINETO GENOME
        return  self.distinct_kmers - sum(self.histo[:coverage_cutoff])


    def set_cutoff(self, cutoff):
        self.kmer_cutoff =  cutoff


    def get_shared_count(self):
        return self.shared_count

    ####DELETE#####
    def percentage_of_kmers(self):
        count = 1
        sum = 0.0
        for v in self.histo_total:
            sum += v
            print count, sum,  "{:.1f}%".format(sum / self.total_kmers * 100.0)
            count += 1

            if count > self.coverage:
                return

    def generate_coverage_histogram(self):
        pass


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
        return (self.name, _arr)

    def compare_to_ardb(self, jfdb):
        pass
