__author__ = 'keith simmon'

import sys, os
import os.path


def clean_string(s):
    text_to_scrub = ["complete sequence",  "whole genome shotgun sequence",
                      "sp.", "strain", ",", "complete genome", "DNA", "chromosome" ]
    if "plasmid" in s:
        s = s[s.find("plasmid") + 8 : ]
    for i in text_to_scrub:
        s = s.replace(i, "")
    return s.replace("   ", " ").replace("  ", " ").strip()



#f you need to modify the taxonomy, edits can be made to the names.dmp and nodes.dmp files in this directory; the gi_taxid_nucl.dmp file will also need to be updated appropriately.

plasmid_path= "/home/ksimmon/reference/strian_typing_resources/bacteria_plasmids/"
files = [f for f in os.listdir(plasmid_path) if f.endswith(".fa")]
plasmid_acc = []
print len(files)
for f in files:
    header = open(plasmid_path + f, "r").readline()
    #print header
    plasmid_acc.append(( header.split("|")[1], header.split("|")[3]))






#root == 1
# 1       |       1       |       no rank |               |       8       |       0       |       1       |       0       |       0       |       0       |       0       |  0|               |

    #Genome  2
        #genome_1
        #genome_2
    #Plasmid 3
        #plasmid_1
        #plasmid_2


enterococcus_path = "/home/ksimmon/reference/strian_typing_resources/enterococcus_genomes/"
files = os.listdir(enterococcus_path)
enterococcus = []
for f in files:
    header = open(enterococcus_path + f, "r").readline()
    #print header.split(" ")[1:3]
    #plasmid_acc.append((header.split("|")[1],header.split("|")[3]))
    gi = header.split("|")[1]
    enterococcus.append((gi, " ".join(header.split(" ")[1:3])),)

acinetobacter_path = "/home/ksimmon/reference/strian_typing_resources/acinetobacter_genomes/"
files = os.listdir(acinetobacter_path)
acinetobacter = []
for f in files:
    header = open(acinetobacter_path + f, "r").readline()
    #print header.split(" ")[1:3]
    #plasmid_acc.append((header.split("|")[1],header.split("|")[3]))
    gi = header.split("|")[1]
    acinetobacter.append((gi,   " ".join(header.split(" ")[1:3])),)

staph_path = "/home/ksimmon/reference/strian_typing_resources/staphyloccus_genomes/"
files = os.listdir(staph_path)
staph = []
for f in files:
    header = open(staph_path + f, "r").readline()
    #print header.split(" ")[1:3],
    #plasmid_acc.append((header.split("|")[1],header.split("|")[3]))
    gi = header.split("|")[1]
    staph.append((gi," ".join(header.split(" ")[1:3])),)





names_str = []
nodes_str = []
gi_taxid_str = []
buffer="\t|\t8\t|\t0\t|\t1\t|\t0\t|\t0\t|\t0\t|\t0\t|\t0\t|"
##create names and nodes
nodes_str.append("{0}\t|\t{1}\t|\t{2}\t|{3}\n".format(1, 1, "no rank", buffer))
nodes_str.append("{0}\t|\t{1}\t|\t{2}\t|{3}\n".format(2, 1, "phylum", buffer))
nodes_str.append("{0}\t|\t{1}\t|\t{2}\t|{3}\n".format(3, 1, "genus", buffer))
nodes_str.append("{0}\t|\t{1}\t|\t{2}\t|{3}\n".format(4, 1, "genus", buffer))
nodes_str.append("{0}\t|\t{1}\t|\t{2}\t|{3}\n".format(5, 1, "genus", buffer))
nodes_str.append("{0}\t|\t{1}\t|\t{2}\t|{3}\n".format(6, 1, "phylum", buffer))


names_str.append("{0}\t|\t{1}\t|\t\t|\t{2}\t|\n".format(1, "all",   "scientific name"))
names_str.append("{0}\t|\t{1}\t|\t\t|\t{2}\t|\n".format(2, "genome",  "scientific name"))
names_str.append("{0}\t|\t{1}\t|\t\t|\t{2}\t|\n".format(3, "Enterococcus_genome",  "scientific name"))
names_str.append("{0}\t|\t{1}\t|\t\t|\t{2}\t|\n".format(4, "Staphylococcus_genome",  "scientific name"))
names_str.append("{0}\t|\t{1}\t|\t\t|\t{2}\t|\n".format(5, "Acinetobacter_genome",  "scientific name"))
names_str.append("{0}\t|\t{1}\t|\t\t|\t{2}\t|\n".format(6, "plasmid", "scientific name"))






starting_taxid = 7


for i, ref in enumerate(enterococcus + acinetobacter + staph):
    taxid, gi, description = i + starting_taxid, ref[0], ref[-1]
    #print taxid, gi, description
    if "Enterococcus" in description:
        nodes_str.append("{0}\t|\t{1}\t|\t{2}\t|{3}\n".format(taxid, 3,"species",  buffer))
        names_str.append("{0}\t|\t{1}\t|\t\t|\t{2}\t|\n".format(taxid, description,   "scientific name"))
        gi_taxid_str.append("{0}\t{1}\n".format(gi, taxid))

    elif "Staphylococcus" in description:
        nodes_str.append("{0}\t|\t{1}\t|\t{2}\t|{3}\n".format(taxid, 4,"species",  buffer))
        names_str.append("{0}\t|\t{1}\t|\t\t|\t{2}\t|\n".format(taxid, description,   "scientific name"))
        gi_taxid_str.append("{0}\t{1}\n".format(gi, taxid))

    elif "Acinetobacter" in description:
        nodes_str.append("{0}\t|\t{1}\t|\t{2}\t|{3}\n".format(taxid, 5,"species",  buffer))
        names_str.append("{0}\t|\t{1}\t|\t\t|\t{2}\t|\n".format(taxid, description,   "scientific name"))
        gi_taxid_str.append("{0}\t{1}\n".format(gi, taxid))

    else:
        print taxid, gi, description

starting_taxid = taxid + 1
for i, ref in enumerate(plasmid_acc):
    taxid, gi, description = i + starting_taxid, ref[0], ref[-1]
    print taxid, gi, description

    nodes_str.append("{0}\t|\t{1}\t|\t{2}\t|{3}\n".format(taxid, 6,"species",  buffer))
    names_str.append("{0}\t|\t{1}\t|\t\t|\t{2}\t|\n".format(taxid, description,   "scientific name"))
    gi_taxid_str.append("{0}\t{1}\n".format(gi, taxid))




base = "/home/ksimmon/reference/strian_typing_resources/kraken_strain_typing/"
wf = open(base + "taxonomy/gi_taxid_nucl.dmp", "w" )
for i in gi_taxid_str:
    wf.write(i)
wf.close()

wf = open(base + "taxonomy/names.dmp", "w" )
for i in names_str:
    wf.write(i)
wf.close()

wf = open(base + "taxonomy/nodes.dmp", "w" )
for i in nodes_str:
    wf.write(i)
wf.close()