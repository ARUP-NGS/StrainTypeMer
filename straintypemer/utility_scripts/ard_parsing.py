from Bio import SeqIO
import sys

_p = "/home/ksimmon/reference/ard/"
sys.stderr.write("Retrieving antibiotic resistance genes\n")
_l = set([])

descriptions = {}
for i in open(_p + "categories.txt"):
    v = i.strip().split("\t")
    name = ".".join(v[0].split(".")[:-1])
    descriptions.update({name : v})

aro_tags = {}
for i in open(_p + "AROtags.txt"):
    v = i.strip().split("\t")
    #print v
    aro_tags.update({v[2]:v[1]})


print aro_tags
for s in SeqIO.parse(_p + "ARmeta-genes.fa", "fasta"):

    id =  s.description.split(" ")[0]
    species = s.description[s.description.rfind("[") + 1:s.description.rfind("]")]
    aro_tag = [i.split(" ")[0] for i in s.description.split(". ") if "ARO:" in i and "ARO:1000001" not in i]
    print id, species, descriptions[id][1], ",".join([aro_tags[tag] for tag in aro_tag])




