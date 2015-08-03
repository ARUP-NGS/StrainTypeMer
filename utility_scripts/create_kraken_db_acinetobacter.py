__author__ = 'keith simmon'

def clean_string(s):
    text_to_scrub = ["complete sequence",  "whole genome shotgun sequence",
                      "sp.", "strain", ",", "complete genome", "DNA", "chromosome" ]
    if "plasmid" in s:
        s = s[s.find("plasmid") + 8 : ]

    for i in text_to_scrub:
        s = s.replace(i, "")

    return s.replace("   ", " ").replace("  ", " ").strip()



#f you need to modify the taxonomy, edits can be made to the names.dmp and nodes.dmp files in this directory; the gi_taxid_nucl.dmp file will also need to be updated appropriately.

base = "/home/ksimmon/reference/strian_typing_resources/kraken_acineto_db/"
references = []
for line in open("/home/ksimmon/reference/strian_typing_resources/kraken_acineto_db/reference_list.txt"):
    references.append( line.strip().split("|") )


#root == 1
# 1       |       1       |       no rank |               |       8       |       0       |       1       |       0       |       0       |       0       |       0       |  0|               |

    #Acinetobacter_genome  2
        #genome_1
        #genome_2
    #Acinetobacter_plasmid 3
        #plasmid_1
        #plasmid_2

names_str = []
nodes_str = []
gi_taxid_str = []
buffer="\t|\t8\t|\t0\t|\t1\t|\t0\t|\t0\t|\t0\t|\t0\t|\t0\t|"
##create names and nodes
nodes_str.append("{0}\t|\t{1}\t|\t{2}\t|{3}\n".format(1, 1, "no rank", buffer))
nodes_str.append("{0}\t|\t{1}\t|\t{2}\t|{3}\n".format(2, 1, "genus", buffer))
nodes_str.append("{0}\t|\t{1}\t|\t{2}\t|{3}\n".format(3, 1, "genus", buffer))

names_str.append("{0}\t|\t{1}\t|\t\t|\t{2}\t|\n".format(1, "all",   "scientific name"))
names_str.append("{0}\t|\t{1}\t|\t\t|\t{2}\t|\n".format(2, "Acinetobacter genome",  "scientific name"))
names_str.append("{0}\t|\t{1}\t|\t\t|\t{2}\t|\n".format(3, "Acinetobacter plasmid", "scientific name"))


for i, ref in enumerate(references):
    taxid, gi, description = i + 4, ref[1], ref[-1]

    if "plasmid" in description:
        name = clean_string(description)
        nodes_str.append("{0}\t|\t{1}\t|\t{2}\t|{3}\n".format(taxid, 3,"species",  buffer))
        names_str.append("{0}\t|\t{1}\t|\t\t|\t{2}\t|\n".format(taxid, name,   "scientific name"))
        gi_taxid_str.append("{0}\t{1}\n".format(gi, taxid))

    elif "genome" in description:
        name = clean_string(description)
        nodes_str.append("{0}\t|\t{1}\t|\t{2}\t|{3}\n".format(taxid, 2,"species", buffer))
        names_str.append("{0}\t|\t{1}\t|\t\t|\t{2}\t|\n".format(taxid, name,   "scientific name"))
        gi_taxid_str.append("{0}\t{1}\n".format(gi, taxid))


# wf = open(base + "taxonomy/gi_taxid_nucl.dmp", "w" )
# for i in gi_taxid_str:
#     wf.write(i)
# wf.close()
#
# wf = open(base + "taxonomy/names.dmp", "w" )
# for i in names_str:
#     wf.write(i)
# wf.close()
#
# wf = open(base + "taxonomy/nodes.dmp", "w" )
# for i in nodes_str:
#     wf.write(i)
# wf.close()