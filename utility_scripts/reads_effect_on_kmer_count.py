import sys
import subprocess
import argparse
from collections import OrderedDict
import numpy as np

"""
Assuming 150 bp per read
Staph/ENTEROCOCCUS == 3M / 150 == 20,000 reads 1X coverage
Acinetobacter 4M/ 150 == 26,000 reads 1X coverage
"""

def get_estimate_coverage(histo):
        """Estimates the coverage for kmers count > 3 times"""
        kmer_count = 0
        #print histo
        try:
            kmer_count =  np.sum([float(v) for k, v in histo.iteritems() if k >= 3])
            if kmer_count == 0 or \
                            float(np.sum(histo.values()[3:])) == 0:
                _cov = 1
            else:
                _cov = np.sum([float(v) * k for k, v in histo.iteritems()] / kmer_count)
        except:
            _cov = 1
        #print kmer_count, histo
        return _cov, kmer_count


parser = argparse.ArgumentParser(description="calculates kmers generated when incrementing read count")

parser.add_argument("-fq", "--fastq_file", help="the fastq file", type=str, required=True)
parser.add_argument("-i", "--increment_step", help="number of reads to chunk at a time", type=int, required=True)
parser.add_argument("-k", "--kmer_length", help="kmer size to analyze", type=int, default=31)
args = parser.parse_args()

fastq_file = args.fastq_file
increment = args.increment_step * 4
kmer_length = str(args.kmer_length)

#get line count
num_lines = sum(1 for line in open(fastq_file))
num_sequences = num_lines / 4


starting_coverage = 5
for num_of_lines in range(increment, num_lines, increment):
    #head
    _histo = OrderedDict()
    cmd = "head -n {0} {1} > temp.fq".format(num_of_lines, fastq_file)

    subprocess.Popen(cmd, shell=True,)
    op = subprocess.Popen(
        ["jellyfish", "count", "-C", "-s", "2G", "-m", kmer_length, "-t", "8", "-o", "count.jf", "temp.fq"],)

    out, err = op.communicate()


    op = subprocess.Popen(["jellyfish", "histo", "count.jf"], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
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

    cov, kmer_count = get_estimate_coverage(_histo)
    print "{0}\t{1}\t{2}\t{3}".format(starting_coverage, int(kmer_count), num_of_lines / 4, "{:.0f}".format(cov))
    starting_coverage += 5


cmd = "rm temp.fq count.jf"
subprocess.Popen(cmd, shell=True)







