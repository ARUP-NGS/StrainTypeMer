from collections import namedtuple
import sys


results_file = "/Users/ksimmon/Box Sync/ARUP/aci_coverage.txt"

results = {}

Values = namedtuple('Values', ['cov', 'kmer_count', 'reads'])
strain = ""
for line in open(results_file):
    #print line
    if ".fq" in line:
        strain = line.split(".")[0]
        results.update({strain : []})
    else:
        line = line.strip().split("\t")
        l = Values(line[0], line[1], line[2])
        results[strain].append(l)

    print line[2]


print results
sys.exit(1)

print ",{0}".format(",".join(sorted(results)))
for i in range(5,200,5):
    print i,
    _r = []
    for strain in sorted(results):
        _r.append([strain] + [(strain, v.kmer_count, v.reads, v. cov) for v in results[strain] if v.cov == str(i)])

    for i in _r:
        if len(i) == 1:
            print ",",
        else:
            print "," + i[1][1],
    print