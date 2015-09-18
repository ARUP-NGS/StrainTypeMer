#!/home/ksimmon/anaconda/bin/python

__author__ = 'keith simmon'

"""
This script will take as input jf databases from strains and perform a comparison

The input of strains should be provided should be provided as a comma separated list.
assume sample name is before underscore
"""
import jf_object as jfobj
import numpy as np
from numpy import sum
import math
import argparse
import os.path
import sys
import subprocess
import collections
from multiprocessing import Process, Queue
import string
import random
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.backends.backend_pdf import PdfPages


## CONFIGURATION INFORMATION ###########################################################################################
jellyfish_path = "/home/ksimmon/bin/jellyfish-2.2.0/bin/jellyfish"
jfobj.jellyfish_path = jellyfish_path
jfobj.ardb_dir = "/home/ksimmon/reference/strian_typing_resources/ARDB/grouped_fastas/jf_files/dump/"
jfobj.ardb_info = "/home/ksimmon/reference/strian_typing_resources/ARDB/class2info.tab"
########################################################################################################################

def chuck_list(_OrderedDict, chunk_size=5):
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
    for i in range(1, 151):
        histo.update({i:0})
    for k, v in  jf_obj.histo.iteritems():
        if k in histo:
            histo[k] += v
        elif k > 150:
            histo[150] += v
    return histo.keys(), histo.values()

def produce_histograms(jf_objects, ylim=350000, xlim=80):
    """
    This will take the list of jellyfish objects (i.e. strains) and create a histogram of the kmers counts

    :param strain_names: [list of strain names]
    :param jf_objects: [list of jf objects]
    :param xlim: the highest value on the x axis
    :param ylim: the highest value on the y axis
    :return: None writes out multiple pdfs with histograms for each strain
    """
    jf_chunks = chuck_list(jf_objects)
    pdf = PdfPages('histograms.pdf')
    for chunk in jf_chunks:
        fig, ax = plt.subplots(len(chunk), sharex=True, sharey=True,)
        fig.set_size_inches(11, 8.5)
        idx = 0
        plt.xlim(2,xlim)
        plt.ylim(0,ylim)

        fig.subplots_adjust(hspace = .3, wspace=.001)

        for name, strain in chunk:
            freq, count = generate_histo_values(strain)
            width = 0.8 # the width of the bars
            ax[idx].bar(freq, count, width, color='g', linewidth=.5)
            ax[idx].set_title(name, fontsize=13)
            ax[idx].spines['top'].set_visible(False)
            ax[idx].spines['right'].set_visible(False)
            ax[idx].yaxis.set_ticks_position('left')
            ax[idx].xaxis.set_ticks_position('bottom')

            #cutoff
            ax[idx].plot([strain.kmer_cutoff, strain.kmer_cutoff], [0, ylim], "--", color='#C24641' ) ##add threshold line
            ax[idx].annotate('mers counted <{0} excluded'.format(strain.kmer_cutoff),
                             xy=(strain.kmer_cutoff, ylim * .95),
                             xytext=(strain.kmer_cutoff + .4, ylim * .95), fontsize=7, color='#C24641', rotation=90)
            ax[idx].add_patch(Rectangle((0, 0), strain.kmer_cutoff, ylim, alpha=.5, facecolor='#C24641', linewidth=0))

            #estimated coverage
            ax[idx].plot([strain.coverage, strain.coverage], [0, ylim], "k--") ##add threshold line
            ax[idx].annotate('{:.1f}x estimated coverage'.format(strain.coverage), xy=(strain.coverage, ylim * .90),
                             xytext=(strain.coverage + .8, ylim * .90), fontsize=8, color='#483C32')

            #estimate genome size
            ax[idx].annotate('Estimated genome size\n{:,} bp'.format(strain.estimate_genome_size(strain.kmer_cutoff)),
                             xy=(xlim * .90, ylim * .80),
                             xytext=(xlim * .90 , ylim * .80), fontsize=8, color='#483C32', horizontalalignment='right')
            idx += 1
            #~~~~~~~~

        ax[len(chunk) / 2].set_ylabel("kmers with a given count (x100)", fontsize=12)
        ax[len(chunk) - 1].set_xlabel("kmer frequency", fontsize=12)

        _xticks = np.arange(0, xlim + 1, 5)
        _xticks[0] = 2
        plt.xticks(_xticks, fontsize=9)
        plt.yticks(np.arange(0, ylim + 1, 50000), np.arange(0, ylim + 1, 50000) / 100, )

        pdf.savefig(transparent=True)
        plt.close()
    pdf.close()
    return


def output_ardb_information(strain_objs, ardb_results):
    ardb_genes_found = {}
    for strain, ardb_result in ardb_results.iteritems():
        for ardb_gene, kmer_count in ardb_result.iteritems():
            #but all records in ardb_found
            if ardb_gene not in ardb_genes_found:
                for strain in ardb_results.iterkeys():
                    ardb_genes_found.update({ardb_gene : collections.OrderedDict({ k : 0 for k in ardb_results.iterkeys() })})
            ardb_genes_found[ardb_gene][strain] = kmer_count

    ardb_info = strain_objs[strain].ardb_info
    plot_count = 1
    pdf = PdfPages('ardb_information.pdf')
    for ardb_gene, results in ardb_genes_found.iteritems():
        fig, ax = plt.subplots(1, sharex=True, sharey=True,)
        fig.set_size_inches(11, 8.5)
        idx = 0
        plt.ylim(0,ardb_info[ardb_gene].max_length)
        fig.subplots_adjust(hspace = 3.5, wspace=.05, bottom=.45)

        N = len(strain_objs)
        width = 0.4
        ind = np.arange(N)
        ax.bar([i for i in range(len(results))] , [i for i in results.itervalues()] , width,
               color='g', linewidth=.5, alpha=.5)

        ax.set_title("{1}".format(ardb_info[ardb_gene].name, ardb_info[ardb_gene].info,
                                  ardb_info[ardb_gene].num_of_sequences), fontsize=10)

        fig.suptitle("{0} (n={2} genes in ARDB)".format(ardb_info[ardb_gene].name, ardb_info[ardb_gene].info,
                                          ardb_info[ardb_gene].num_of_sequences), fontsize=15, y=.95)

        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.yaxis.set_ticks_position('left')
        ax.xaxis.set_ticks_position('bottom')

        ax.plot([0, len(strain_objs)], [ardb_info[ardb_gene].min_length, ardb_info[ardb_gene].min_length],
                "--", color='#C24641' ) ##add threshold line
        ax.annotate('Smallest gene size in ARDB',
                          xy=(.3, ardb_info[ardb_gene].min_length),
                          xytext=(.3, ardb_info[ardb_gene].min_length * 1.1), fontsize=9, color='#C24641',)

        ax.plot([0, len(strain_objs)], [ardb_info[ardb_gene].avg_length, ardb_info[ardb_gene].avg_length],
                "--", color='b' ) ##add threshold line
        ax.annotate('Average gene size in ARDB',
                          xy=(3.3, ardb_info[ardb_gene].avg_length),
                          xytext=(3.3, ardb_info[ardb_gene].avg_length * .95), fontsize=9, color='b',)


        ax.plot([0, len(strain_objs)], [ardb_info[ardb_gene].median_length, ardb_info[ardb_gene].median_length],
                "--", color='#C24641' ) ##add threshold line
        ax.annotate('Median gene size in ARDB',
                          xy=(.3, ardb_info[ardb_gene].median_length),
                          xytext=(.3, ardb_info[ardb_gene].median_length * .95), fontsize=9, color='#C24641',)

        ax.plot([0, len(strain_objs)], [ardb_info[ardb_gene].max_length, ardb_info[ardb_gene].max_length],
                "--", color='#C24641' ) ##add threshold line
        ax.annotate('Max gene size in ARDB',
                          xy=(.3, ardb_info[ardb_gene].max_length),
                          xytext=(.3, ardb_info[ardb_gene].max_length * .95), fontsize=9, color='#C24641',)

        plt.xticks(ind+(width/2), [label[0:25] for label in strain_objs.keys()], rotation=270)
        pdf.savefig(transparent=True)
        plt.close()
    pdf.close()
    return



def count_kmers(q, merged_jf_obj, this_jf_object, cutoff):
    """
    This function facilitates the queuing of the jellyfish counting
    Takes a jellyfish object merged from all the strains and calls the gets kmer count for the queried strain
    :param q: queue
    :param merged_jf_obj: merged count of all strains
    :param this_jf_object: strain jf object
    :return: None (puts in queue)    """
    q.put(this_jf_object.get_kmer_count(merged_jf_obj, cutoff))
    return


def main():
    """
    The main method which takes a list jf counts (strains) compares them.
    :return: status
    """
    parser = argparse.ArgumentParser(description="the script take multiple jellyfish counts and compares the strains")

    parser.add_argument("-c", "--cutoff", help="The number of kmers to analyze [Default = None]", type=int,
              default=None)

    parser.add_argument("--coverage_cutoff", help="percent of genome coverage to set kmer filters [DEFAULT .20 if "  + \
                                                  "coverage is 30 [(30 * .20) = 6] kmers with a count < 5 will " + \
                                                  "be ignored for cooresponding strain", type=float, default=.20)

    parser.add_argument("-t", "--cpus",
              help="The number of cpus to use when counting kmers in strains [Default is the len of the strain list]",
              type=int, default=2)

    parser.add_argument( "--no_kmer_filtering", help="Do not filter kmers based on coverage",
              action="store_true", default=False)

    parser.add_argument("-k", "--kmer_reference",
              help="instead of merging strains use kmer reference set for comparison",
              type=str, default=None)

    parser.add_argument("--do_not_output_histograms",
                        help="This will prevent the output of the PDF files containing the histograms",
                        default=True, action="store_false")

    parser.add_argument('jf_files', nargs='+', help='jellyfish files for each strain')

    args = parser.parse_args()
    no_kmer_filtering = args.no_kmer_filtering
    cutoff = args.cutoff
    cpus = args.cpus
    coverage_cutoff = args.coverage_cutoff
    jf_files = args.jf_files
    kmer_reference = args.kmer_reference
    output_histogram = args.do_not_output_histograms

    compare_strains(jf_files=jf_files, no_kmer_filtering=no_kmer_filtering, cutoff=cutoff, cpus=cpus,
                    coverage_cutoff=coverage_cutoff,kmer_reference=kmer_reference, output_histogram=output_histogram)
    return

def compare_strains(jf_files, no_kmer_filtering, cutoff, cpus, coverage_cutoff, kmer_reference, output_histogram):
    ####################################################################################################################
    # CHECK ARGUMENTS AND ATTRIBUTES
    ####################################################################################################################
    if kmer_reference is not None and os.path.isfile(kmer_reference) is False:
        sys.stderr.write("kmer reference file does not exist: {0}\n".format(kmer_reference))
        sys.exit(11)

    NUM_OF_STRAINS = len(jf_files) #calculate the number of samples
    if cpus == None or cpus > NUM_OF_STRAINS:
        cpus = NUM_OF_STRAINS


    file_paths = []
    strain_objs = collections.OrderedDict()  ##KEEP RECORD OF FILE PATHS and CREATE THE STRAIN OBJS FOR EACH STRAIN
    for f in jf_files:
        if os.path.exists(f) == False: ##MAKE SURE THE FILE EXISTS
            sys.stderr.write("file does not exist: {0}\n".format(f))
            sys.exit(1)
        file_paths.append(f)
        strain_name = "".join(os.path.basename(f).split(".")[0])

        strain_objs.update({strain_name : jfobj.jf_object(strain_name, f)})

    if len(strain_objs) != NUM_OF_STRAINS: ##MAKE SURE THE NAMES ARE UNIQUE
        sys.stderr.write("strain names are not unique:\n{0}\n".format(", ".join(strain_objs.keys())))
        sys.exit(2)

    ####################################################################################################################
    temp_file =  "/tmp/tmp_{0}".format(''.join(random.choice(string.ascii_uppercase) for i in range(8)))
    jf_temp_file = temp_file + ".jf"
    dump_temp_file = temp_file + ".txt"
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #  COMPARE STRAINS DIRECTLY
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if kmer_reference is None:
        try:
            subprocess.check_call([jellyfish_path, "merge", '-L', "2", "-o", jf_temp_file] + file_paths )
        except:
            sys.stderr.write("Error in running jellyfish merge\n")
            sys.exit(3)
        # dump kmers out to a file [Technically I should be able to remove this however a jellyfish bug causes the loss
        # of canonicalization for merged sets ] dumping out to txt file avoids this.
        try:
            subprocess.check_call([jellyfish_path, "dump", '-c',  "-t", "-o", dump_temp_file, jf_temp_file] )
        except:
            sys.stderr.write("Error in running jellyfish merge\n")
            sys.exit(4)

        os.remove(jf_temp_file)
        merged_jf = dump_temp_file
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #  USE USER SUPPLIED KMER LIST  [kmer\tcount]
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    else:
        merged_jf = kmer_reference



    ########################################################################################################################
    #START THE WORK
    ########################################################################################################################
    ardb_results = {}
    count_table = [[] for i in range(len(strain_objs)) ]
    counter = 0
    attach = count_table.append

    #~~~~~~~~~~~~~~~ MULTIPROCESS THIS MOFO ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    q = Queue() #TODO EXPLORE PIPES
    #create the job queue
    jobs = strain_objs.values()
    num_of_strains_counted = 0
    current_processes = []


    for cpu in range(cpus):
        _obj = jobs.pop()
        p = Process(target=count_kmers, args=(q, merged_jf, _obj, cutoff,), name=_obj.name)
        current_processes.append(p)

    #start the jobs for the correct number of cpus
    while len(current_processes) != 0:
        current_processes.pop().start()

    #keep the queue moving
    while len(jobs) != 0:
        _name, _arr, ardb = q.get() #PAUSES UNTIL A JOBS RETURN  #AM I WIATING FOR THE FIRST OBJECT
        ardb_results.update({_name : ardb})
        num_of_strains_counted += 1

        sys.stderr.write("{0}\tof\t{1}\tstrains processed\n".format(num_of_strains_counted, len(strain_objs)))


        count_table[strain_objs.keys().index(_name)]=np.array(_arr)
        #start next job
        p = Process(target=count_kmers, args=( q, merged_jf, jobs.pop(), cutoff, ),name=_obj.name)
        p.start()
    #nothing else to start

    #wait until the queue returns 'ALL THE THINGS'
    while num_of_strains_counted != len(strain_objs): ##finished processing
        _name, _arr, ardb = q.get() #PAUSES UNTIL A JOBS RETURN
        ardb_results.update({_name : ardb})
        num_of_strains_counted += 1

        sys.stderr.write("{0}\tof\t{1}\tstrains processed\n".format(num_of_strains_counted, len(strain_objs)))
        count_table[strain_objs.keys().index(_name)]=np.array(_arr)
    q.close()
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    #cleanup temp files
    if os.path.isfile(dump_temp_file):
        os.remove(dump_temp_file)


    sys.stdout.write(" COVERAGE INFORMATION ".center(80, "-") + "\n")
    sys.stdout.write("calculations only include kmers counted > 3 times\n")


    #determine kmer cutoff for each strain ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    cutoff_table = [] #this will hold the cutoff
    if no_kmer_filtering:
        for k, v in strain_objs.iteritems():
            sys.stdout.write("Strain: {:s}\t Coverage Estimate: {:.1f}\n".format(k, v.coverage))
            #cutoff_table.append(0)
            v.set_cutoff(0)


    else:
        for k, v in strain_objs.iteritems():
            sys.stdout.write("Strain: {:s}\t Coverage Estimate: {:.1f}\n".format(k, v.coverage))
            cutoff_table.append(math.ceil(v.coverage * coverage_cutoff))
            if math.ceil(v.coverage * coverage_cutoff) <= 3:
                sys.stderr.write(" WARNING ".center(80, "-") + "\n")
                sys.stderr.write("Strain: {0} has low coverage\n".format(k))
                sys.stderr.write("Calculated cutoff is {0}\n".format(int(math.ceil(v.coverage * coverage_cutoff))))
                if (math.ceil(v.coverage * coverage_cutoff)) < 3:
                    sys.stderr.write("Changing kmer cutoff to 3\n")
                    v.set_cutoff(3)
                sys.stderr.write("If estimated genome size is lower than expected consider repeating\n")
                sys.stderr.write("".center(80, "-") + "\n\n")
                v.set_cutoff(int(math.ceil(v.coverage * coverage_cutoff)))
            else:
                v.set_cutoff(int(math.ceil(v.coverage * coverage_cutoff)))
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    sys.stdout.write(" STRAIN STATS ".center(80, "-") + "\n")
    for k, v in strain_objs.iteritems():
        sys.stdout.write("Strain: {:s}\n".format(k))
        sys.stdout.write("\tInferred genome size: {:.0f}  [filtering kmers counted <= {:.0f} times]\n".format(
                          v.estimate_genome_size(v.kmer_cutoff), v.kmer_cutoff ))
    sys.stdout.write("".center(80, "-") + "\n\n")


    #FILTER TABLE CREATED
    cutoff_cov = np.array( [v.kmer_cutoff for v in strain_objs.values()] )

    ##FILTER COUNT TABLE
    count_table_filtered = []
    total_kmers = 0.0
    kept_kmers = 0.0

    for kmer_count in np.array(count_table).T:
        #print kmer_count
        total_kmers += 1
        if (kmer_count - cutoff_cov).max() >= 0:
            count_table_filtered.append(kmer_count.clip(0,1))
            kept_kmers += 1

    sys.stdout.write(" FILTER INFO ".center(80, "-") + "\n")
    sys.stdout.write("{:.0f}\tkmers counted (n > 1) in all strains\n".format(total_kmers))
    sys.stdout.write("{:.0f}\tkmers included in analysis\n".format(kept_kmers))
    sys.stdout.write("".center(80, "-") + "\n\n")

    count_table_filtered = np.array(count_table_filtered).T #column become row and rows become columns.
    strain_keys = strain_objs.keys()


    # WRITE OUT HISTOGRAMS
    if output_histogram:
        produce_histograms(strain_objs)

    #DETERMINE RELATIONSHIPS~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    similarity_dict = collections.OrderedDict()
    for i in range(len(count_table_filtered)):
        similarity_dict.update({strain_keys[i] : collections.OrderedDict()})
        sum_1 = np.sum(count_table_filtered[i])  #sum of kmers strain 1
        for j in range(i+1, len(count_table_filtered)):
             sum_2 = np.sum(count_table_filtered[j]) #sum of kmers strain 2
             intersection = sum( (count_table_filtered[i] + count_table_filtered[j]) == 2 )
             total_kmers =  sum( (count_table_filtered[i] + count_table_filtered[j]) > 0 )
             smallest = sum_1
             if sum_2 < smallest:
                 smallest = sum_2
             #print sum_1, sum_2, total_kmers, intersection, smallest
             #print intersection, total_kmers
             similarity_dict[strain_keys[i]].update({
                                                     strain_keys[j] :
                                                         (float(intersection) / total_kmers * 100,
                                                          float(intersection) / total_kmers * 100,
                                                          total_kmers,
                                                          smallest,
                                                          )
                                                    })

    #PRINT SIMILARITY TABlE~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    sys.stdout.write("[SIMILARITY TABLE]\n")
    delimeter = ","
    _str = delimeter + delimeter.join(strain_keys) + "\n"
    for i in range(len(strain_keys)):
        _str += strain_keys[i] + delimeter
        for j in range(len(strain_keys)):
            if i == j:
                _str +=  "100" + delimeter
            elif i < j:
                _str += "{:.1f}{:s}".format(similarity_dict[strain_keys[i]][strain_keys[j]][0], delimeter)
            elif i > j:
                _str += "{:.1f}{:s}".format(similarity_dict[strain_keys[j]][strain_keys[i]][1], delimeter)
        _str  = _str[:-1] + "\n"
    sys.stdout.write(_str)
    sys.stdout.write("[SIMILARITY TABLE END]\n")
    #ADD DENOMINATOR TO OUTPUT
    sys.stdout.write("\n\n")
    sys.stdout.write("[DENOMINATOR TABLE]\n")
    _str = delimeter + delimeter.join(strain_keys) + "\n"
    for i in range(len(strain_keys)):
        _str += strain_keys[i] + delimeter
        for j in range(len(strain_keys)):
            if i == j:
                _str += "{0}{1}".format(np.sum(count_table_filtered[j]), delimeter)
            elif i < j:
                _str += "{:.1f}{:s}".format(similarity_dict[strain_keys[i]][strain_keys[j]][2], delimeter)
            elif i > j:
                _str += "{:.1f}{:s}".format(similarity_dict[strain_keys[j]][strain_keys[i]][3], delimeter)
        _str  = _str[:-1] + "\n"
    sys.stdout.write(_str)
    sys.stdout.write("[DENOMINATOR TABLE END]\n")
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    output_ardb_information(strain_objs, ardb_results)
    return
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

if __name__ == "__main__":
    main()