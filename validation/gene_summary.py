table_file = "/Users/331-SimmonkLPTP/Documents/data/strain_typing/acinetobacter_baumannii/info/ACI_ARD.txt"


start_analysis = False

results = {}
for line in open(table_file):
   #print(line.strip())

    if "-- STRAIN STATS ------" in line:
        start_analysis = True

    if start_analysis:
        if line.startswith("Strain:"):
            strain_name = line.split(":")[-1].strip()
            results.update({strain_name : []})
        else:

            if "ARD GENE:" in line:
                _r = []
                line = line.split(":")
                gene_tag = line[2].replace("Covered", "").strip()
                covered = float(line[3].strip().split("%")[0])
                _r.append(gene_tag)
                #_r.append(covered)
            if "count mean" in line:
                coverage = line.split(";")[-1].split(" ")[1]
                #_r.append(coverage)
            if "Description:" in line:
                desc = line.split(":")[-1].strip()
                #_r.append(desc)
                results[strain_name].append(", ".join(_r))

d = sorted(results, key=results.get)
for k in d:
    print(k, sorted(results[k]))
    # for i in results[k]:
    #     print(i[0])
