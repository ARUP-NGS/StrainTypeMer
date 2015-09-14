__author__ = 'ksimmon'

import argparse

parser = argparse.ArgumentParser(description="the script take multiple jellyfish counts and compares the strains")
parser.add_argument('vcf_files', nargs='+', help='jellyfish files for each strain', type=argparse.FileType("r"))

args = parser.parse_args()
#no_kmer_filtering = args.no_kmer_filtering
#cutoff = args.cutoff
#cpus = args.cpus
vcf_files = args.vcf_files

vcf_results = []
coverage_stats = []
for f in vcf_files:
    _results = {}
    _cov = 0
    for line in f:
        if line[0] != "#":
            line = line.strip().split("\t")
            pos, alt, score = int(line[1]), line[4], float(line[5])
            if score > 200:
                info = { i.split("=")[0] : i.split("=")[1] for i in line[7].split(";") if "=" in i}
                coverage = int(info["DP"])
                no_of_ref_reads, no_of_alt_reads = 0, 0
                try:
                    no_of_ref_reads, no_of_alt_reads =  [int(v) for v in info["DPR"].split(",")] ##variant
                except:
                    try:
                        no_of_ref_reads = int(info["DPR"])  ## non_variant
                    except: ## multic alleles
                        pass
                        # try:
                        # #     no_of_ref_reads, no_of_alt_reads_1, no_of_alt_reads_2 =  [int(v) for v in info["DPR"].split(",")] ##variant
                        # #     if no_of_alt_reads_1 > no_of_alt_reads_2:
                        # #         no_of_alt_reads = no_of_alt_reads_1
                        # #         alt = alt.split(",")[0]
                        # #     else:
                        # #         no_of_alt_reads = no_of_alt_reads_2
                        # #         alt = alt.split(",")[1]
                        # except:
                        #     no_of_ref_reads, no_of_alt_reads = 0, 0 ## mask this position
                        #     print line
                mapping_qual= info["MQ"]
                try:
                    mapping_qual = int(mapping_qual)
                except:
                    mapping_qual = 0

                #passing filter
                if mapping_qual > 10 and (no_of_alt_reads >= 15 or no_of_ref_reads >= 15):
                    if no_of_alt_reads > no_of_ref_reads:
                        if float( no_of_alt_reads ) / (no_of_alt_reads + no_of_ref_reads) > .80:
                            if pos not in _results:
                                _results.update({pos : alt})
                                _cov += no_of_alt_reads + no_of_ref_reads
                            else:
                                _results.pop(pos)
                    else:
                        if float( no_of_ref_reads ) / (no_of_alt_reads + no_of_ref_reads) > .80:
                            if pos not in _results:
                                _results.update({pos : "."})
                                _cov += no_of_alt_reads + no_of_ref_reads
                            else:
                               _results.pop(pos)
    #print "strain counted"
    print "{0}\t{1}\tavg read depth".format(f.name.split("/")[-1].split("_")[0], "{:.1f}".format(float(_cov)/ len(_results)))
    vcf_results.append(_results)



#intersection_set = set(vcf_results[0]).intersection(set(vcf_results[1]))


#tot, dif, = 0.0, 0.0


results_diffs =        [[0 for i in range(len(vcf_files))] for j in range(len(vcf_files))]
results_intersection = [[0 for i in range(len(vcf_files))] for x in range(len(vcf_files))]


for i in range(len(vcf_files)):
    for j in range(len(vcf_files)):

        if i==j:
            results_diffs[i][j] = 0
            results_intersection[i][j] = (len(vcf_results[i]), len(vcf_results[i]), len(vcf_results[i]),)


        else:
            tot, dif, = 0.0, 0.0
            #print( i,j)
            intersection_set = set(vcf_results[i]).intersection(set(vcf_results[j]))
            for pos in intersection_set:
                tot += 1
                if vcf_results[i][pos] != vcf_results[j][pos]:
                    #print pos, vcf_results[i][pos], vcf_results[j][pos]
                    dif += 1


            results_diffs[i][j] = results_diffs[j][i] = int(dif)
            results_intersection[i][j] = ( len(intersection_set), len(vcf_results[i]), len(vcf_results[j]), )
            results_intersection[j][i] = ( len(intersection_set), len(vcf_results[i]), len(vcf_results[j]), )

            #print "~" * 80
            #print "Sample_1:\t{0}\nSample_2:\t{1}".format(vcf_files[i].name.split("/")[-1].split("_")[0], vcf_files[j].name.split("/")[-1].split("_")[0])
            #print "s1:\t{0}\tpositions passing filters\ns2:\t{1}\tpositions passing filters".format(len(vcf_results[i]), len(vcf_results[j]))
            #print "{0}\tintersecting positions passing filters".format(len(intersection_set))
            #print "{0}\tgenomic differences".format(int(dif))
            #print "{0}%\tpercent of genomic differences".format("{:.2f}".format(dif/tot * 100))
            #print "~" * 80
            #print

_str =","
for f in vcf_files:
    _str += f.name.split("/")[-1].split("_")[0] + ","

_str = _str[:-1] + "\n"
for i, f in enumerate(vcf_files):
    _str += f.name.split("/")[-1].split("_")[0] + ","
    for j,v in enumerate(results_diffs[i]):
        _str += str(v) + ","
    _str = _str[:-1] + "\n"

print
print
print _str
print
print

_str = ","
for f in vcf_files:
    _str += f.name.split("/")[-1].split("_")[0] + ","
_str = _str[:-1] + "\n"
for i, f in enumerate(vcf_files):
    _str += f.name.split("/")[-1].split("_")[0] + ","
    for j,v in enumerate(results_intersection[i]):
        _str += str("i:{0};s1:{1};s2:{2}".format( v[0],v[1],v[2], )) + ","
    _str = _str[:-1] + "\n"
print _str