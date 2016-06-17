import os, sys
import matplotlib.pyplot as plt
import numpy as np

sample = 'mrsa'
sample = 'enterococcus'
#sample = 'acinetobacter'
ylimit = 4000000

plt.figure(figsize=(11, 8))
ax = plt.subplot(111)
ax.spines["top"].set_visible(False)
ax.spines["bottom"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.spines["left"].set_visible(False)
ax.get_xaxis().tick_bottom()
ax.get_yaxis().tick_left()
plt.ylim(0, ylimit)
plt.xlim(0, 2500000)

plt.yticks(range(0, ylimit + 100000, 500000), [str(float(x/1000000)) + "e6" for x in
                                               range(0, ylimit + 100000, 500000)], fontsize=14)
plt.xticks(range(0, 2600000, 500000), ["{0:,.0f}".format(x) for x in range(0, 2600000, 500000)], fontsize=14)
plt.xticks(fontsize=14)



plt.tick_params(axis="both", which="both", bottom="off", top="off",
                labelbottom="on", left="off", right="off", labelleft="on")


base_path="/Users/331-SimmonkLPTP/Box Sync/ARUP/rarefaction/{0}/".format(sample)
files_to_plot = [os.path.join(base_path, i) for i in os.listdir(base_path) if ".txt" in i]

print(files_to_plot)

#total_kmers, distinct_kmers, read_count, coverage, kmer_filter
all_results = {}
for f in files_to_plot:
    x , y = [0], [0]
    for line in open(f, "r"):
        line = line.strip().split("\t")
        x.append(int(line[2]))
        y.append(int(line[1]))
    s = os.path.basename(f).split("_")[0][:-1]
    replicate = os.path.basename(f).split("_")[0]
    distinct_kmers, coverage = int(line[1]), float(line[3])

    if s in all_results:
        all_results[s].append((replicate, distinct_kmers, coverage,))
    else:
        all_results.update({s: [(replicate, distinct_kmers, coverage,)]})

    plt.plot(x, y, '-', alpha=.5)

for k,v in all_results.items():
    try:
        avg, dif = np.average([v[0][1], v[1][1]]), abs(v[0][1] - v[1][1])
        var = "{:.2f}%".format(float(dif) / avg * 100)
        print("{0}\t{1}\t{2}\t{3}\t{4}\t{5}".format(k, var, int(avg), dif, v[0][2], v[1][2]))
    except:
        pass

for y in range(0, ylimit + 100000, 500000):
    plt.plot([2500000, 0], [y, y], "--", lw=0.5, color="black", alpha=0.3)

plt.ylabel("Number of distinct kmers", fontsize=16)
plt.xlabel("Number of reads", fontsize=16)
plt.title("Distinct kmers acquired per reads sampled")
plt.savefig(os.path.join(base_path, "{0}_nextseq_rarefaction.pdf".format(sample)), bbox_inches="tight", dpi=800, format='pdf')

