import numpy as np
from collections import OrderedDict
table_file = "/Users/331-SimmonkLPTP/Documents/data/strain_typing/acinetobacter_baumannii/info/mrsa_tables.csv"

def get_stats(a):
    return np.mean(a), np.median(a), np.max(a), np.min(a), len(a)



tables = OrderedDict()
kmer_table = None
relationship_table = {}
is_first = True
for line in open(table_file):
    #print(line.strip())
    if "TABLE:" in line:
        if kmer_table is not None:
            tables.update({table_name : [relationship_table , kmer_table]})
            is_first = False
        kmer_table = {}
        #relationship_table = {}
        table_name = line.split(",")[0].split(":")[-1].strip()
    elif line[0] == ",": #column labels
        col_labels = line.strip().split(",")
    else: # DATA
        _row = line.strip().split(",")
        #print(_row)
        is_pfge = True
        for i, value in enumerate(_row):
            if i == 0:
                row_label = value
            elif value == "":
                is_pfge = False
            elif is_pfge:
                if is_first:
                    name = "{0}:{1}".format(row_label, col_labels[i])
                    if name in relationship_table:
                        relationship_table[name][0] = value
                    else:
                        relationship_table.update({name : value})
            else:
                name = "{0}:{1}".format(col_labels[i], row_label,)
                if name in kmer_table:
                    kmer_table[name][1] = float(value)
                else:
                    kmer_table.update({name: float(value)})
else:
    tables.update({table_name: [relationship_table , kmer_table]})

all = []
for name, table in tables.items():
    results = {name : {"C": [], "P": [], "U": [], "I" : []}}

    for comparison, value in table[0].items():
        results[name][value].append(table[1][comparison])

    all.append(results)


for comparison in all:

    for name, results in comparison.items():

        identical = np.array(results["I"])
        close = np.array(results["C"])
        prop = np.array(results["P"])
        unrelated = np.array(results["U"])

        new_array = np.array([results["I"], results["C"],results["P"],results["U"]])

        longest = 0
        for i in new_array:
            if len(i) > longest:
                longest = len(i)

        print(longest, len(new_array[0]))
        print(name)
        print("I\tC\tP\tU")

        for i in range(longest):
            I, C, P, U = "", "", "", ""

            try:
                I = new_array[0][i]
            except:
                pass
            try:
                C = new_array[1][i]
            except:
                pass
            try:
                P = new_array[2][i]
            except:
                pass
            try:
                U = new_array[3][i]
            except:
                pass

            print("{0}\t{1}\t{2}\t{3}".format(I,C,P,U))


        print()




        i = get_stats(identical)
        c = get_stats(close)
        p = get_stats(prop)
        u = get_stats(unrelated)

        # print(name)
        # print("category,mean,median,max,min,n")
        # print("I,{0}".format(",".join(["{:.1f}".format(j) for j in i])))
        # print("C,{0}".format(",".join(["{:.1f}".format(j) for j in c])))
        # print("P,{0}".format(",".join(["{:.1f}".format(j) for j in p])))
        # print("U,{0}".format(",".join(["{:.1f}".format(j) for j in u])))
        # print()
