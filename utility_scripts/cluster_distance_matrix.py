__author__ = 'ksimmon'

import scipy
import pylab
import scipy.cluster.hierarchy as sch
import argparse
import numpy as np

parser = argparse.ArgumentParser(description="the script take multiple jellyfish counts and compares the strains")
parser.add_argument("-m", "--matrix", help="the distance matrix", type=argparse.FileType('r'), required=True)


args = parser.parse_args()
matrix = args.matrix

x_labels = [ i.split("_")[0] for i in matrix.next().split(",")[1:]]
y_labels = []
data = []
for line in matrix:
    line = line.strip().split(",")
    y_labels.append(line[0].split("_")[0])
    data.append([float(i) for i in line[1:]])

print data

D = np.array(data)

# Compute and plot first dendrogram.
fig = pylab.figure(figsize=(8,8))
ax1 = fig.add_axes([0.07,0.1,0.15,0.6], frame_on=False)
Y = sch.linkage(D, method='single')
Z1 = sch.dendrogram(Y, orientation='right', labels=y_labels, truncate_mode='lastp', )
ax1.set_xticks([])
ax1.set_yticks([])

# Compute and plot second dendrogram.
ax2 = fig.add_axes([0.3,0.75,0.6,0.2], frame_on=False)
Y = sch.linkage(D, method='single', )
print Y
Z2 = sch.dendrogram(Y, labels=x_labels, color_threshold=None, count_sort=False, )
ax2.set_xticks([])
ax2.set_yticks([])

# Plot distance matrix.
axmatrix = fig.add_axes([0.3,0.1,0.6,0.6], frame_on=False)
idx1 = Z1['leaves']
idx2 = Z2['leaves']
D = D[idx1,:]
D = D[:,idx2]
im = axmatrix.matshow(D, aspect='auto', origin='lower', cmap=pylab.cm.Greens, )
print x_labels
axmatrix.set_xticks(range(len(x_labels)))
axmatrix.set_xticklabels([x_labels[i] for i in idx1], minor=False)
axmatrix.xaxis.set_label_position("top")
axmatrix.xaxis.tick_top()

pylab.xticks(rotation=-90, fontsize=8)

axmatrix.set_yticks(range(len(x_labels)))
axmatrix.set_yticklabels([y_labels[i] for i in idx2], minor=False)
axmatrix.yaxis.set_label_position('left')
axmatrix.yaxis.tick_left()

pylab.yticks(fontsize=8)

axcolor = fig.add_axes([0.94,0.1,0.02,0.6],frame_on=False)

# Plot colorbar.
pylab.colorbar(im, cax=axcolor)

fig.show()
fig.savefig('/home/ksimmon/dendrogram_cookson.png')