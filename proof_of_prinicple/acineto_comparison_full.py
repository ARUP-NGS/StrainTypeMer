__author__ = 'ksimmon'


import numpy as np
import math
#import brewer2mpl
import pandas as pd
import scipy
import pylab
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import scipy.spatial.distance as distance
import scipy.cluster.hierarchy as sch

base= "/home/ksimmon/data/strain_typing/Acineto_combined/jellyfish_31/dump/recount/MERGED/"
merged_file="acineto_merged_31bp.fa"

order = []
NUMBER_OF_STRAINS = 21

strain_values = {}
total_kmers = {}

COMPLETELY_SHARED = 0

_key = { 1: "A1",
          2 : "A2",
         4: "A3",
         8: "A4_2",
         16: "A4",
         32: "A5_2",
         64: "A5",
         128 :"A6",
         256 : "A7",
         512: "A8",
         1024 : "A9",
         2048 : "AC1",
         4096 : "AC2",
         8192 : "AC3",
         16384 : "AC4",
         32768 : "AC5",
         65536 : "AC6",
         131072 : "AC7",
         262144 : "AC8",
         524288 : "AC9_2",
         1048576 : "AC9" }


for i in range(NUMBER_OF_STRAINS):
    offset = int(math.pow(2,i))
    strain_values.update({offset : {} })
    total_kmers.update({offset : 0 })
    order.append(offset)
    COMPLETELY_SHARED+=offset
    for j in range(i+1,NUMBER_OF_STRAINS):
        offset_inner = int(math.pow(2,j))
        strain_values[offset].update({offset_inner:0})

wf = open('/home/ksimmon/Documents/atrib_1000_remove_shared.txt', "w")

count = 0
for line in open(base + merged_file):
    #kmer = line.strip().split("\t")[0]
    flag = int(line.strip().split("\t")[1])
    if flag != COMPLETELY_SHARED:
        _s = ""
        for i in range(len(order)):
            if order[i] & flag:
                _s += "1,"
            else:
                _s += "0,"
        wf.write(_s[:-1] + "\n")

        count += 1
        if count == 10000:
            break



print COMPLETELY_SHARED
count = 0
for line in open(base + merged_file):
    #kmer = line.strip().split("\t")[0]
    flag = int(line.strip().split("\t")[1])
    if flag != COMPLETELY_SHARED:

        for i in range(len(order)):
            if order[i] & flag:

                total_kmers[order[i]] += 1
                for j in range(i+1,len(order)):
                    if order[j] & flag:
                        strain_values[order[i]][order[j]] += 1
        count += 1
        if count == 10000:
            break
################print out formatted results###########################

distance_matrix = []

print "".center(5) + ",",
for i in order[:]:
    print _key[i].center(5) + ",",
print

distance_matrix = []
for i in order:
    _arr = []
    print _key[i].center(5) + "," ,
    for j in order[:]:
        #print i, j, strain_values[i][j], total_kmers[i], total_kmers[j]
        if j in strain_values[i]:
            _arr.append(  float(strain_values[i][j] + strain_values[i][j]) / float(total_kmers[i] + total_kmers[j]))
            print "{:.1f},".format((float(strain_values[i][j] + strain_values[i][j])) / float(total_kmers[i] + total_kmers[j]) * 100).center(5),
        elif i == j:
            print "100".center(5) + ",",
            _arr.append(1)
        else:
            _arr.append( float(strain_values[j][i] + strain_values[j][i]) / float(total_kmers[i] + total_kmers[j]))
            print "{:.1f},".format((float(strain_values[j][i] + strain_values[j][i])) / float(total_kmers[i] + total_kmers[j]) * 100).center(5),
    distance_matrix.append(_arr)
    print ""

distance_matrix = np.array(distance_matrix)



#######################MAKE THIS A STANDALONE MODULE##############################################3


_lab = []
for i in order:
    _lab.append(_key[i])

print _lab

# _labels_1 = []
# for i in idx1:
#     _labels_1.append(_lab[i])
#
# print
# _labels_2 = []
# for i in idx2:
#     _labels_2.append(_lab[i])
#pairwise_dists_rows = distance.squareform(distance.pdist(distance_matrix))
#print "number of rows: {0}".format(distance_matrix.shape[0])
#print "size of matrix: {0}".format(pairwise_dists_rows.shape)

#pairwise_dists_cols = distance.squareform(distance.pdist(distance_matrix.T))

#clusters_rows = sch.linkage(pairwise_dists_rows, method='average') #UPMGA
clusters_rows = sch.linkage(distance_matrix, method='average') #UPMGA
clusters_cols = sch.linkage(distance_matrix, method='single') #NEAREST POINT





fig = plt.figure()
heatmapGS = gridspec.GridSpec(2,2,wspace=0.0,hspace=0.0,width_ratios=[0.25,1], height_ratios=[0.25,1])

#col DENDRO
col_denAX = fig.add_subplot(heatmapGS[0,1])
col_denD = sch.dendrogram(clusters_cols, color_threshold=np.inf)
col_ax = col_denAX.get_axes()
col_ax.set_xticks([])
col_ax.set_yticks([])

##ROW DENDRO
row_denAX = fig.add_subplot(heatmapGS[1,0])
row_denD = sch.dendrogram(clusters_rows, color_threshold=np.inf, orientation="right")
row_ax = row_denAX.get_axes()
row_ax.set_xticks([])
row_ax.set_yticks([])

#HEATMAP

heatmapAX = fig.add_subplot(heatmapGS[1,1])
idx1 = col_denD['leaves']
idx2 = row_denD['leaves']

D = distance_matrix[idx1,:]
D = D[:,idx2]

axi = heatmapAX.imshow(D,interpolation='nearest',aspect='auto', origin='lower',  cmap=plt.cm.Blues)
ax = axi.get_axes()
ax.set_xticks([])
ax.set_yticks([])



###ADD LABELS
##row
heatmapAX.yaxis.set_ticks_position("right")
idx1_labels = [_lab[i]  for i in idx1]
heatmapAX.set_yticks(np.arange(distance_matrix.shape[0]))
heatmapAX.set_yticklabels(idx1_labels)


##col
heatmapAX.xaxis.set_ticks_position("bottom")
idx2_labels = [_lab[i]  for i in idx2]
heatmapAX.set_xticks(np.arange(distance_matrix.shape[0]))
xls = heatmapAX.set_xticklabels(idx2_labels)
##rotate
for label in xls:
    label.set_rotation(90)

for l in heatmapAX.get_xticklines() + heatmapAX.get_yticklines():
    l.set_markersize(0)


###COLOR BAR###
scale_cbGSSS = gridspec.GridSpecFromSubplotSpec(1,2,subplot_spec=heatmapGS[0,0], wspace=0.0, hspace=0.0)

#axcolor = fig.add_axes([0.91,0.1,0.02,0.6])
axcolor = fig.add_subplot(scale_cbGSSS[0,0])
cb = fig.colorbar(axi,axcolor)
cb.set_label("percent id")
cb.ax.yaxis.set_ticks_position('left')
cb.ax.yaxis.set_label_position('left')
start, end = cb.ax.get_ylim()
print start, end

#cb.ax.yaxis.set_ticks(np.arange(start + .2, end - .2, .2))
plt.tick_params(labelsize=8)
cb.outline.set_linewidth(0)

fig.tight_layout()
fig.savefig('/home/ksimmon/Desktop/dendro_10000_remove_shared.pdf')
plt.show()

