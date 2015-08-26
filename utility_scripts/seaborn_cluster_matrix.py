__author__ = 'ksimmon'
import pandas as pd
import seaborn as sns
sns.set(font="monospace")

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




# Load the brain networks example dataset
df = sns.load_dataset("brain_networks", header=[0, 1, 2], index_col=0)


print df

# Select a subset of the networks
used_networks = [1, 5, 6, 7, 8, 11, 12, 13, 16, 17]
used_columns = (df.columns.get_level_values("network")
                          .astype(int)
                          .isin(used_networks))
df = df.loc[:, used_columns]

# Create a custom palette to identify the networks
network_pal = sns.cubehelix_palette(len(used_networks),
                                    light=.9, dark=.1, reverse=True,
                                    start=1, rot=-2)
network_lut = dict(zip(map(str, used_networks), network_pal))

# Convert the palette to vectors that will be drawn on the side of the matrix
networks = df.columns.get_level_values("network")
network_colors = pd.Series(networks).map(network_lut)

# Create a custom colormap for the heatmap values
cmap = sns.diverging_palette(h_neg=210, h_pos=350, s=90, l=30, as_cmap=True)

# Draw the full plot
sns.clustermap(df.corr(), row_colors=network_colors, linewidths=.5,
               col_colors=network_colors, figsize=(13, 13), cmap=cmap)
sns.plt.show()