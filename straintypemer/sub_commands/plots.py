import matplotlib
matplotlib.use("Agg")
import numpy as np
import scipy.cluster.hierarchy as sch
import collections
import matplotlib.pyplot as pylab
from matplotlib.patches import Rectangle
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import colors
matplotlib.rcParams['lines.linewidth'] = 1

def adjust_font_size(matrix_length):
    """
    Returns a int with a font size that will fit in a matrix of X length
    :param matrix_length: the size one edge of the matrix
    :return: (int) the suggested font size for given matrix size
    """
    adjusted_size = 0
    if matrix_length < 7:
        adjusted_size = 25
    if matrix_length >= 7:
        adjusted_size = 20
    if matrix_length >= 9:
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
        adjusted_size = 1
    if matrix_length > 100:
        adjusted_size = 0
    return adjusted_size

def generage_matrix(x_labels, y_labels, data, output_prefix, kmer_count, vmin=50,
                    identical=None, possibly_related=None, different=None):
    """
    This function saves a PDF file containing a clustered matrix (nearest neighbor).

    :param x_labels: The labels for the x axis
    :param y_labels: The labels for the y axis (most of the time same as the x_labels)
    :param data: The matrix (for kmer typing this should be expressed as percent identity (0-100))
    :param output_prefix: Prefix for the file_name
    :param kmer_counts: the number of kmers analyzed {strain_label : num_kmers}
    :param vmin: default is 50, this will automatically be set to the lowest value for
    :return: None (writes file out)
    """
    matplotlib.rcParams['lines.linewidth'] = .8
    D = np.array(data, dtype= 'float')

    if identical is not None:
        cmap = colors.ListedColormap(['#992600', '#b3b300', '#006600'])
        bounds = [0, possibly_related, identical, 100]
        norm = colors.BoundaryNorm(bounds, cmap.N)
    else:
        cmap = pylab.cm.RdYlGn
        norm = None

    #truncate labels
    fig = pylab.figure()
    fig.set_size_inches(11, 8.5)

    # Compute and plot first dendrogram. [LEFT]
    #ax1 = fig.add_axes([0.058,0.1,0.155,0.6], frame_on=False, )
    ax1 = fig.add_axes([0.045,0.1,0.155,0.6], frame_on=False, )
    Y = sch.linkage(D, method='weighted')
    Z1 = sch.dendrogram(Y, orientation='left', labels=y_labels, color_threshold=0, )  # color_list=['k'] )
    ax1.set_xticks([])
    ax1.set_yticks([])

    # Compute and plot second dendrogram. [TOP]
    #ax2 = fig.add_axes([0.26,0.815,0.6,0.1], frame_on=False)
    ax2 = fig.add_axes([0.26,0.775,0.6,0.1], frame_on=False,)
    Y = sch.linkage(D, method='weighted')
    Z2 = sch.dendrogram(Y, labels=x_labels, count_sort=False, color_threshold=0, no_plot=True) # color_list=['k'])
    ax2.set_xticks([])
    ax2.set_yticks([])

    # Plot distance matrix.
    axmatrix = fig.add_axes([0.26,0.1,0.6,0.6], frame_on=False)
    idx1 = Z1['leaves'] #returns clustered ordered labels
    idx2 = Z2['leaves']

    D = D[idx1,:]
    D = D[:,idx2]

    #see if the minumum is too high
    if np.min(D) < vmin:
        vmin = np.min(D)


    # minor hacking to create the minor ticks, this is need to overlay the grid
    # probably a little bit of over kill but I think it better shows each grids as a distinct observation
    matplotlib.rcParams['lines.linewidth'] = 0.8
    locs = np.arange(len(x_labels))
    adjusted_locs = []
    for loc in locs:
        adjusted_locs.append(loc + .5)

    for axis in [axmatrix.xaxis, axmatrix.yaxis]:
        axis.set_ticks(adjusted_locs, minor=True, )
        axis.set(ticks=locs + 0.5, ticklabels=x_labels, )

    im = axmatrix.matshow(D, aspect='auto', origin='lower', cmap=cmap, norm=norm,
                          interpolation='nearest', vmin=vmin, vmax=100, )
    matplotlib.rcParams['lines.linewidth'] = 1.5
    #modifiy x ticks and labels
    axmatrix.set_xticks(range(len(x_labels)))
    axmatrix.set_xticklabels(["{0:s}".format((x_labels[i].replace(".fq", "").replace(".fastq", "").replace(".gz", "")[:10] + "\n" +
                                              x_labels[i].replace(".fq", "").replace(".fastq", "").replace(".gz", "")[10:20]).replace("_"," "),
                                              kmer_count[x_labels[i]]) for i in idx2], minor=False,
                                              multialignment='center', linespacing=1)


    axmatrix.xaxis.set_label_position("top")
    if len(x_labels) > 100:
        pylab.xticks(rotation=-90, fontsize=2, )
    else:
        pylab.xticks(rotation=-90, fontsize=8,)
        pylab.yticks(fontsize=3,)

    #make sure minor ticks are one
    axmatrix.tick_params(axis='both', which='minor', right='on', top='on', bottom='on', color="w", width=.02)

    #turn off the major ticks
    axmatrix.tick_params(axis='both', which='major', right='off', top='off', bottom='off', left='off', width=.02)

    #modify the y tick and labels

    axmatrix.set_yticks(range(len(y_labels)))
    axmatrix.set_yticklabels(["{:s}".format((y_labels[i].replace(".fq", "").replace(".fastq", "").replace(".gz", "")[:10]
                                             + "\n" +
                                             y_labels[i].replace(".fq", "").replace(".fastq", "").replace(".gz", "")[10:20]).replace("_"," "),
                                             kmer_count[y_labels[i]]) for i in idx1], minor=False,
                                             multialignment='center', linespacing=1)

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

    if len(x_labels) > 100:
        pylab.yticks(fontsize=2)
    else:
        fs = int(adjust_font_size(len(data)) / 2)
        pylab.yticks(fontsize=fs)
        pylab.xticks(fontsize=adjust_font_size(len(data)))




    #~~~ Add annotations to each cell
    font_size = adjust_font_size(len(data))

    if len(x_labels) < 100:
        for x in range(len(x_labels)):
            for y in range(len(y_labels)):
                #modifiy formatting to make sure value fits in square
                val = "{:.1f}".format(D[x][y])
                _color = 'k'
                if "100" in val:
                    val = '100'
                    _color = '#C3C300'

                #use white font color if value is greater than 99

                # if float(D[x][y]) >= possibly_related:
                #     _color = 'black'
                # else:
                #     _color = "w"
                #
                #
                # if float(D[x][y]) >= identical:
                #     _color = 'w'


                #annotate this mother
                axmatrix.annotate(val, xy=(x, y), horizontalalignment='center', verticalalignment='center',
                                  fontsize=font_size, color=_color )


    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #color scale bar
    if identical is None:
        axcolor = fig.add_axes([0.88,0.1,0.02,0.6], frame_on=False)
        pylab.colorbar(im, cax=axcolor)

    #make a grid on the minor axis
    axmatrix.grid(True, which='minor', linestyle='-', color="w", linewidth=0.1)

    #filename
    if output_prefix != "":
        pdf = "{0}_{1}".format(output_prefix,'matrix.png',)
    else:
        pdf = 'matrix.png'

    fig.savefig(pdf, format="png", bbox_inches="tight", dpi=400)
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
        _l.append(list(_OrderedDict.items())[i: i + chunk_size])
    return _l

def generate_histo_values(jf_obj):
    """
    returns the freq and count of 300 kmers

    :param jf_obj: will retrieve the histo object
    :return: frequency count [1-150] and count at each frequency
    """
    histo = collections.OrderedDict()
    for i in range(1, 501):
        histo.update({i: 0})
    for k, v in  jf_obj.histo.items():
        if k in histo:
            histo[k] += v
        elif k > 500:
            histo[500] += v
    return list(histo.keys()), list(histo.values()), max(list(histo.values())[10: -2]) * 1.5

def produce_histograms(jf_objects, output_prefix):
    """
    This will take the list of jellyfish objects (i.e. strains) and create a histogram of the kmers counts

    :param strain_names: [list of strain names]
    :param jf_objects: [list of jf objects]
    :return: None writes out multiple pdfs with histograms for each strain
    """
    jf_chunks = chunk_list(jf_objects, )
    if output_prefix != "":
        pdf = PdfPages("{0}_{1}".format(output_prefix,'histograms.png'))
    else:
        pdf = PdfPages('histograms.png')

    ylim = 0
    #plt.use("Agg")
    for chunk in jf_chunks:

        fig = pylab.figure()
        fig.set_size_inches(11, 8.5)

        idx = 0

        pylab.subplots_adjust(hspace = .4, wspace=.001)
        this_plot = 1
        for name, strain in chunk:
            ax = fig.add_subplot(5 , 1, this_plot)

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
                        xytext=(xlim * .90, ylim * .60), fontsize=8, color='#483C32', horizontalalignment='right')

        #     if this#~~~~~~~~
            pylab.xticks(fontsize=9)
            pylab.yticks(fontsize=9)
            if this_plot == 4:
                ax.set_ylabel("kmers with a given count", fontsize=12)

        if this_plot < 5:
            ax.set_ylabel("kmers with a given count", fontsize=12)

        ax.set_xlabel("kmer frequency", fontsize=12)
        pdf.savefig(transparent=True)
        pylab.close()
    pdf.close()
    return
