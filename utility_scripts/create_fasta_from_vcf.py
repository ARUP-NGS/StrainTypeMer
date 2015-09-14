__author__ = 'ksimmon'

import argparse
import sys

parser = argparse.ArgumentParser(description="the script take multiple jellyfish counts and compares the strains")
parser.add_argument("-v", '--vcf_file', help='jellyfish files for each strain', type=argparse.FileType("r"))

args = parser.parse_args()
#no_kmer_filtering = args.no_kmer_filtering
#cutoff = args.cutoff
#cpus = args.cpus
vcf_file = args.vcf_file

count = 0
vcf_results = []
coverage_stats = []
_results = {}
_cov = 0
sequence = {}
for i in vcf_file:

    if i[0] != "#":
        count += 1
        line = i.strip().split("\t")
        pos = line[1]
        ref = line[3]
        alt = line[4]
        score = line[5]
        info = { i.split("=")[0] : i.split("=")[1] for i in line[7].split(";") if "=" in i}
        coverage = int(info["DP"])
        no_of_ref_reads, no_of_alt_reads = 0, 0
        try:
            no_of_ref_reads, no_of_alt_reads =  [int(v) for v in info["DPR"].split(",")] ##variant
        except:
            try:
                no_of_ref_reads = int(info["DPR"])  ## non_variant
            except: ## multiaellic alleles
                pass
        mapping_qual= info["MQ"]
        try:
            mapping_qual = int(mapping_qual)
        except:
            mapping_qual = 0

        #print pos, ref, alt, coverage
        #passing filter

        if mapping_qual > 10 and (no_of_alt_reads >= 10 or no_of_ref_reads >= 10):
            if no_of_alt_reads > no_of_ref_reads:
                if float( no_of_alt_reads ) / (no_of_alt_reads + no_of_ref_reads) > .80:
                    if pos not in _results: ### ALT
                        _results.update({pos : alt})
                        _cov += no_of_alt_reads + no_of_ref_reads
                    else:
                        sys.stderr.write(_results[pos])
                        print line
                        print
                        _results.pop(pos)
            else:
                if float( no_of_ref_reads ) / (no_of_alt_reads + no_of_ref_reads) > .80:
                    if pos not in _results:
                        _results.update({pos : "."})
                        _cov += no_of_alt_reads + no_of_ref_reads
                    else:
                        _results.pop(pos)
                        print pos

        else:
            _results.update({pos : "X"})
            #print pos, ref, alt, coverage, no_of_alt_reads, no_of_ref_reads



count = 1
for k, v in _results.iteritems():
    if k != count:
        print k, count
    count +=1