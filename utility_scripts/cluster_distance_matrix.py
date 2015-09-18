__author__ = 'ksimmon'

#import scipy
import pylab
import scipy.cluster.hierarchy as sch
import argparse
import math
import numpy as np
#import sys

parser = argparse.ArgumentParser(description="the script take multiple jellyfish counts and compares the strains")
parser.add_argument("-m", "--matrix", help="the distance matrix", type=argparse.FileType('r'), required=True)
parser.add_argument("-o", "--output", help="the distance matrix", type=str, required=True)

args = parser.parse_args()
matrix = args.matrix
output = args.output

x_labels = [] #[ i.split("_")[0] for i in matrix.next().strip().split(",")[1:]]
y_labels = []
data = []

is_similarity_table = False
for line in matrix:
    if "[SIMILARITY TABLE]" in line:
        #print line
        is_similarity_table = True
    if "[SIMILARITY TABLE END]" in line:
        #print line
        is_similarity_table = False

    if is_similarity_table and "[SIMILARITY TABLE" not in line:

        if str(line).startswith(","):
            x_labels = [ i.split("_")[0] for i in line.strip().split(",")[1:]]
        else:
            line = line.strip().split(",")
            y_labels.append(line[0].split("_")[0])

            data.append([float(i) for i in line[1:]])




#print x_labels
#print y_labels

D = np.array(data, dtype= 'float')

if np.max(D) > 100:
    pass

    #D = (1 - (D / np.max(D))) * 100

    #D = np.log10(D + 0.01)

# Compute and plot first dendrogram.
fig = pylab.figure()
fig.set_size_inches(11, 8.5)
ax1 = fig.add_axes([0.058,0.1,0.15,0.6], frame_on=False)
Y = sch.linkage(D, method='weighted')
Z1 = sch.dendrogram(Y, orientation='right', labels=y_labels, color_threshold=0, color_list=['k'] )

ax1.set_xticks([])
ax1.set_yticks([])

# Compute and plot second dendrogram.
ax2 = fig.add_axes([0.26,0.75,0.6,0.2], frame_on=False)
Y = sch.linkage(D, method='weighted', )
#print Y
Z2 = sch.dendrogram(Y, labels=x_labels, count_sort=False, color_threshold=0, color_list=['k'])
ax2.set_xticks([])
ax2.set_yticks([])

# Plot distance matrix.
axmatrix = fig.add_axes([0.26,0.1,0.6,0.6], frame_on=False)
idx1 = Z1['leaves']
idx2 = Z2['leaves']
#print idx1

#print D[idx1,:]
D = D[idx1,:]
D = D[:,idx2]
im = axmatrix.matshow(D, aspect='auto', origin='lower', cmap=pylab.cm.RdYlGn ,  )
print im.get_contains()
#print x_labels
axmatrix.set_xticks(range(len(x_labels)))
axmatrix.set_xticklabels([x_labels[i] for i in idx1], minor=False)
axmatrix.xaxis.set_label_position("top")
axmatrix.xaxis.tick_top()
axmatrix.minorticks_off()

axmatrix.tick_params(
    axis='both',
    which='both',      # both major and minor ticks are affected
    left='off',
    right='off',# ticks along the bottom edge are off
    top='off',

        # ticks along the top edge are off
   ) # labels along the bottom edge are off


pylab.xticks(rotation=-90, fontsize=8)

axmatrix.set_yticks(range(len(x_labels)))
axmatrix.set_yticklabels([y_labels[i] for i in idx2], minor=False)
axmatrix.yaxis.set_label_position('left')
axmatrix.yaxis.tick_left()



for x in xrange(len(x_labels)):
    for y in xrange(len(y_labels)):
        val = str(D[x][y])
        if val == '100.0':
            val = '100'

        axmatrix.annotate(val, xy=(y, x),
                    horizontalalignment='center',
                    verticalalignment='center', fontsize=5)

pylab.yticks(fontsize=8)

axcolor = fig.add_axes([0.88,0.1,0.02,0.6],frame_on=False)

# Plot colorbar.
pylab.colorbar(im, cax=axcolor)

#fig.show()
fig.savefig(output)