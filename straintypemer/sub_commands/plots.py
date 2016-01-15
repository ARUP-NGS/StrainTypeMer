import matplotlib
matplotlib.use("Agg")
import numpy as np
import pylab
import scipy.cluster.hierarchy as sch
import collections
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.backends.backend_pdf import PdfPages

def adjust_font_size(matrix_length):
    """
    Returns a int with a font size that will fit in a matrix of X length
    :param matrix_length: the size one edge of the matrix
    :return: (int) the suggested font size for given matrix size
    """
    adjusted_size = 0
    if  matrix_length < 7:
        adjusted_size = 25
    if matrix_length >=7:
        adjusted_size = 20
    if matrix_length >=9:
        adjusted_size = 18
    if matrix_length > 10:
        adjusted_size = 14
    if matrix_length > 12:
        adjusted_size = 12
    if matrix_length > 15:
        adjusted_size = 10
    if matrix_length > 20:
        adjusted_size = 8
    if matrix_length > 25:
        adjusted_size = 6
    if matrix_length > 30:
        adjusted_size = 4
    if matrix_length > 40:
        adjusted_size = 3
    if matrix_length > 50:
        adjusted_size = 0
    return adjusted_size



def generage_matrix(x_labels, y_labels, data, output_prefix, strain_lengths, vmin=50):
    """
    This function saves a PDF file containing a clustered matrix (nearest neighbor).

    :param x_labels: The labels for the x axis
    :param y_labels: The labels for the y axis (most of the time same as the x_labels)
    :param data: The matrix (for kmer typing this should be expressed as percent identity (0-100))
    :param output_prefix: Prefix for the file_name
    :param vmin: default is 50, this will autamatically be set to the lowest value for
    :return: None (writes file out)
    """
    D = np.array(data, dtype= 'float')

    #truncate labels
    #x_labels_trun = [lab[:7] for lab in x_labels]
    #y_labels_trun = [lab[:7] for lab in y_labels]

    fig = pylab.figure()
    fig.set_size_inches(11, 8.5)

    # Compute and plot first dendrogram. [LEFT]
    ax1 = fig.add_axes([0.058,0.1,0.115,0.6], frame_on=False, )
    Y = sch.linkage(D, method='weighted')
    Z1 = sch.dendrogram(Y, orientation='right', labels=y_labels, color_threshold=0, color_list=['k'] )
    ax1.set_xticks([])
    ax1.set_yticks([])

    # Compute and plot second dendrogram. [TOP]
    ax2 = fig.add_axes([0.26,0.815,0.6,0.1], frame_on=False)
    Y = sch.linkage(D, method='weighted', )
    Z2 = sch.dendrogram(Y, labels=x_labels, count_sort=False, color_threshold=0, color_list=['k'])
    ax2.set_xticks([])
    ax2.set_yticks([])

    # Plot distance matrix.
    axmatrix = fig.add_axes([0.26,0.1,0.6,0.6], frame_on=False)
    idx1 = Z1['leaves'] #returns clustered ordered lables
    idx2 = Z2['leaves']

    D = D[idx1,:]
    D = D[:,idx2]

    #see if the minumum is too high
    if np.min(D) < vmin:
        vmin = np.min(D)

    im = axmatrix.matshow(D, aspect='auto', origin='lower', cmap=pylab.cm.RdYlGn,
                         interpolation='nearest',vmin=vmin, vmax=100,)

    # minor hacking to create the minor ticks, this is need to overlay the grid
    # probably a little bit of over kill but I think it better shows each grids as a distinct observartion
    locs = np.arange(len(x_labels))
    for axis in [axmatrix.xaxis, axmatrix.yaxis]:
        axis.set_ticks(locs + 0.5, minor=True)
        axis.set(ticks=locs, ticklabels=x_labels)

    #modifiy x ticks and labels
    axmatrix.set_xticks(range(len(x_labels)))
    axmatrix.set_xticklabels(["{:s}\n(n={:,} kmers)".format(x_labels[i][:20].replace("_"," "),
                                    strain_lengths[x_labels[i]]) for i in idx2 ], minor=False, multialignment='center',
                             linespacing=1)

    axmatrix.xaxis.set_label_position("top")
    pylab.xticks(rotation=-90, fontsize=6,)

    #make sure minor ticks are one
    axmatrix.tick_params(axis='both', which='minor', right='on', top='on', bottom='on', color="w")

    #turn of the major ticks
    axmatrix.tick_params(axis='both', which='major', right='off', top='off', bottom='off', left='off')

    #modify the y tick and labels
    axmatrix.set_yticks(range(len(y_labels)))
    axmatrix.set_yticklabels( ["{:s}\n(n={:,} kmers)".format(y_labels[i][:20].replace("_"," "),
                                    strain_lengths[y_labels[i]]) for i in idx1 ], minor=False, multialignment='center',
                              linespacing=1)

    for i, label in enumerate(axmatrix.xaxis.get_ticklabels()):
        if i%2 == 0:
            label.set_color('#5C6161')
        else:
            label.set_color("#000000")

    for i, label in enumerate(axmatrix.yaxis.get_ticklabels()):
        if i%2 == 0:
            label.set_color('#5C6161')
        else:
            label.set_color("#000000")

    pylab.yticks(fontsize=6)
    #~~~ Add annotations to each cell
    font_size = adjust_font_size(len(data))
    for x in xrange(len(x_labels)):
        for y in xrange(len(y_labels)):
            #modifiy formating to make sure value fits in square
            val = "{:.1f}".format(D[x][y])
            _color = 'k'
            if "100" in val:
                val = '100'
                _color = '#C3C300'

            #use white font color if value is greater than 90

            if float(D[x][y]) > 95.0:
                _color = 'w'

            if float(D[x][y]) > 99.9:
                _color = '#C3C300'
            #annotate this mother
            axmatrix.annotate(val, xy=(x, y), horizontalalignment='center', verticalalignment='center',
                              fontsize=font_size, color=_color )
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #color scale bar
    axcolor = fig.add_axes([0.88,0.1,0.02,0.6], frame_on=False)
    pylab.colorbar(im, cax=axcolor)

    #make a grid on the minor axis
    axmatrix.grid(True, which='minor', linestyle='-', color="w", linewidth=0.5)

    #filename
    if output_prefix != "":
        pdf = "{0}_{1}".format(output_prefix,'matrix.pdf')
    else:
        pdf = 'matrix.pdf'

    fig.savefig(pdf)
    return

def chunk_list(_OrderedDict, chunk_size=5):
    """
    little method that helps to break the number of strains into chunks so that the histograms are printed on a
    reasonably scaled PDF

    :param _list: 1d list
    :param chunk_size: chunks size
    :return: 2d array with equal chunks [last chunk may be smaller]
    """
    _l = []
    for i in range(0, len(_OrderedDict), chunk_size):
        _l.append(_OrderedDict.items()[i : i + chunk_size])
    return _l

def generate_histo_values(jf_obj):
    """
    returns the freq and count of 300 kmers

    :param jf_obj: will retrieve the histo object
    :return: frequency count [1-150] and count at each frequency
    """
    histo = collections.OrderedDict()
    for i in range(1, 351):
        histo.update({i:0})
    for k, v in  jf_obj.histo.iteritems():
        if k in histo:
            histo[k] += v
        elif k > 350:
            histo[350] += v
    return histo.keys(), histo.values(), max(histo.values()[3:-2]) * 1.5

def produce_histograms(jf_objects, output_prefix):
    """
    This will take the list of jellyfish objects (i.e. strains) and create a histogram of the kmers counts

    :param strain_names: [list of strain names]
    :param jf_objects: [list of jf objects]
    :return: None writes out multiple pdfs with histograms for each strain
    """
    jf_chunks = chunk_list(jf_objects, )
    if output_prefix != "":
        pdf = PdfPages("{0}_{1}".format(output_prefix,'histograms.pdf'))
    else:
        pdf = PdfPages('histograms.pdf')

    ylim = 0
    #plt.use("Agg")
    for chunk in jf_chunks:

        fig = plt.figure()
        fig.set_size_inches(11, 8.5)

        idx = 0

        plt.subplots_adjust(hspace = .4, wspace=.001)
        this_plot = 1
        for name, strain in chunk:
            ax = fig.add_subplot(5,1,this_plot)

            freq, count, ylim = generate_histo_values(strain)
            xlim = strain.coverage * 2
            width = 0.8 # the width of the bars

            ax.set_title(name, fontsize=13)
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)
            ax.yaxis.set_ticks_position('left')
            ax.xaxis.set_ticks_position('bottom')

            ax.bar(freq, count, width, color='g', linewidth=.5,)

            this_plot += 1

            ax.set_xlim(2, xlim)
            ax.set_ylim(0, ylim)
            #cutoff
            ax.plot([strain.kmer_cutoff, strain.kmer_cutoff], [0, ylim], "--", color='#C24641' )#add threshold line
            ax.annotate('<{0} excluded'.format(strain.kmer_cutoff),
                              xy=(strain.kmer_cutoff, ylim * .95),
                              xytext=(strain.kmer_cutoff + .4, ylim * .95),
                              fontsize=6, color='#C24641', rotation=90)

            ax.add_patch(Rectangle((0, 0), strain.kmer_cutoff, ylim, alpha=.5, facecolor='#C24641', linewidth=0))
        #
        #     #estimated coverage
            ax.plot([strain.coverage, strain.coverage], [0, ylim], "k--") ##add threshold line
            ax.annotate('{:.1f}x estimated coverage'.format(strain.coverage), xy=(strain.coverage, ylim * .90),
                              xytext=(strain.coverage + .8, ylim * .90), fontsize=8, color='#483C32')
        #
        #     #estimate genome size
            ax.annotate('Estimated genome size\n{:,} bp'.format(strain.estimate_genome_size(strain.kmer_cutoff)),
                             xy=(xlim * .90, ylim * .80),
                             xytext=(xlim * .90 , ylim * .60), fontsize=8, color='#483C32', horizontalalignment='right')
        #
        #     if this#~~~~~~~~
            plt.xticks(fontsize=9)
            plt.yticks(fontsize=9)
            if this_plot == 4:
                ax.set_ylabel("kmers with a given count", fontsize=12)


        if this_plot < 5:
            ax.set_ylabel("kmers with a given count", fontsize=12)

        ax.set_xlabel("kmer frequency", fontsize=12)




        pdf.savefig(transparent=True)
        plt.close()
    pdf.close()
    return
