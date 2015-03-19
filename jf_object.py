__author__ = 'ksimmon'

import jellyfish
import subprocess
import sys
import numpy as np

class jf_object:
    """contains information about the jellyfish db for a single strain"""

    def __init__(self, name, path):
        self.name = name
        self.path = path
        self.qf = jellyfish.QueryMerFile(path)

        self.get_stats()
        self.get_histo()
        self.get_estimate_coverage()
        self.estimate_genome_size()
        self.shared_count = None
        #self.percentage_of_kmers()


    def get_stats(self):
        op = subprocess.Popen(["jellyfish", "stats", self.path], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
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
        op = subprocess.Popen(["jellyfish", "histo", "--threads", "4", "--full", "--high", str(self.max_count - 1),
                              "--low", "2", self.path], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
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
        self.coverage = float(np.sum(self.histo_total[3:])) / float(np.sum(self.histo[3:]))


    def get_histo_plot(self):
        pass
        #TODO

    def estimate_genome_size(self):  ####GET THE FREQUENCY OF KMERS REPEATED IN ACINETO GENOME
        print self.name, self.coverage
        print  self.distinct_kmers
        print  self.distinct_kmers - sum(self.histo[:5])
        print  self.distinct_kmers - sum(self.histo[:10])
        print  self.distinct_kmers - sum(self.histo[:15])

        print

    def kmer_count(self, file_path):
        _arr = []
        attach = _arr.append
        count = 0
        with open(file_path) as f:
            for l in f:
                count +=1
                kmer = l.strip().split("\t")[0]
                mer = jellyfish.MerDNA(kmer)
                mer.canonicalize()
                attach(self.qf[mer])
                if count > 1000:
                    break
        self.shared_count = _arr
        return _arr

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