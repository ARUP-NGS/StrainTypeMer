#!/home/ksimmon/anaconda/bin/python

"""
This application will take as fasta/fastq sequences from strains and perform a comparison in kmer space.
The input of strains should be provided as a positional argument
"""
from straintypemer.sub_commands.compare import compare
from straintypemer.sub_commands.update_mlst_resources import update_mlst_resources
from straintypemer.sub_commands.count import count
from straintypemer.sub_commands.plot_output import plot_output
import argparse

########################################################################################################################
########################################################################################################################


def arguments():
    """
    Parses user supplied arguments for StrainTypeMer sub-commands
    :return: named tuple of arguments
    """
    parser = argparse.ArgumentParser(prog="StrainTypeMer", formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description="the script take multiple jellyfish counts and compares the strains")

    subparsers = parser.add_subparsers(help='sub-command help', dest='subparser_name', title="sub-commands")

    subparsers.add_parser('update_mlst', help='update mlst resources')

    parser_fastq = subparsers.add_parser('compare', help='compare strains')

    parser_fastq.add_argument("-t", "--cpus",
                              help="The number of cpus to use [Default: 1]",
                              type=int, default=1)

    parser_fastq.add_argument("--no_kmer_filtering",
                              help="Do not filter kmers based on coverage; useful when comparing reference sequences",
                              action="store_true", default=False)

    parser_fastq.add_argument("-q", "--qual_score", help="the phred score to filter bases", default=0, type=int)

    parser_fastq.add_argument("-k", "--kmer_reference",
                              help="Use kmer reference set for comparison (e.g. plasmid, core genome, pan genome kmers;"
                                   " [ONLY COMPARE KMERS IN THIS FILE])", type=str, default=None,)

    parser_fastq.add_argument("-r", "--inverse_kmer_reference",
                              help="Use kmer reference set for comparison (e.g. plasmid, core genome, pan genome kmers;"
                                   " [ONLY COMPARE KMERS NOT IN THIS FILE])",
                              type=str, default=None, )

    # hide this in help; allows user to supply same inverse and kmer reference
    parser_fastq.add_argument("-kr", "-rk", help=argparse.SUPPRESS, type=str, default=None)

    parser_fastq.add_argument("--coverage_cutoff",
                              help="percent of genome coverage to set kmer filters [DEFAULT .20 if coverage is 30 " +
                                   "[(30 * .20) = 6] kmers with a count < 5 will be ignored for corresponding strain",
                              type=float, default=0.2)

    parser_fastq.add_argument("--do_not_output_histograms",
                              help="This will prevent the output of the PDF files containing the histograms",
                              default=True, action="store_false")

    parser_fastq.add_argument("--do_not_output_matrix",
                              help="This will prevent the output of the PDF files containing the matrix", default=True,
                              action="store_false")

    parser_fastq.add_argument("--no_pdfs", help="Output will only go to stdout", default=False, action="store_true")

    parser_fastq.add_argument("-o", "--output_prefix", help="appends a prefix to the output files", default="")

    parser_fastq.add_argument("-gz", "--gzipped", help="flag to indicate fastq/a files are gzipped [Default False]",
                              action="store_true", default=False)

    parser_fastq.add_argument("-pwf", "--pairwise_kmer_filter",
                              help="evaluate non-shared kmers in closely related strains",
                              action="store_true", default=False)

    parser_fastq.add_argument("-ard", "--include_ard_comparison", help="include comparison with ard genes",
                              action="store_true", default=False)

    parser_fastq.add_argument("--rapid_mode", help="analyze a subset of kmers", action="store_true",
                              default=False)

    parser_fastq.add_argument('-jf', "--jf_input", help=argparse.SUPPRESS, action="store_true",
                              default=False)

    parser_fastq.add_argument('-kf', "--keep_files", help=argparse.SUPPRESS, action="store_false",
                              default=True)

    parser_fastq.add_argument('fq_files', nargs='+',
                              help='fastq files for each strain (fq1 OR fq1;label OR fq1,fq2;label will be matching '
                                   'string of the two files OR fq1,fq2;label included "[NF]" at the end to prevent '
                                   'kmer filter (useful when adding a reference genome)')

    parser_plot = subparsers.add_parser('plot', help='plot results from standard out')

    parser_plot.add_argument("-i", "--input", help="the input file", type=argparse.FileType("r"), required=True)
    parser_plot.add_argument("-o", "--output_prefix", help="the output prefix", type=argparse.FileType("w"),
                             required=True)

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # parser_count = subparsers.add_parser('count', help='counts and creates a strain_object that contains information'
    #                                                    'about the strain',)
    #
    # parser_count.add_argument("-t", "--cpus",
    #                           help="The number of cpus to use when counting kmers in strains [Default == 2]",
    #                           type=int, default=2)
    #
    # parser_count.add_argument("-gz", "--gzipped", help="flag to indicate fastq_files are gzipped [Default False]",
    #                           action="store_true", default=False)
    #
    # parser_count.add_argument("--coverage_cutoff",
    #                           help="percent of genome coverage to set kmer filters [DEFAULT .15 if coverage is 30 " +
    #                                "[(30 * .15) = 4.5 ~= 5] kmers with a count < 5 will be ignored for corresponding "
    #                                "strain",
    #                           type=float, default=.15)
    #
    #  parser_count.add_argument("--no_kmer_filtering", help="Do not filter kmers based on coverage",
    # action="store_true",
    #                           default=False)
    #
    # parser_count.add_argument("-q", "--qual_score", help="the phred score to filter bases", default=0, type=int)
    #
    # parser_count.add_argument("-l", "--label", help="the label to attach to the file", default=None)
    #
    # parser_count.add_argument("-o", "--out", help="the output file")
    #
    # parser_count.add_argument('fq_files', nargs='+',
    #                           help='fastq files for a strain')

    args = parser.parse_args()
    return args
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


def main():
    args = arguments()
    if args.subparser_name == "count":
        count(
            fq_files=args.fq_files, gzipped=args.gzipped,
            no_kmer_filtering=args.no_kmer_filtering,
            coverage_cutoff=args.coverage_cutoff,
            qual_filter=args.qual_score, cpus=args.cpus,
            label=args.label, out=args.out
        )

    if args.subparser_name == "compare":
        if args.no_pdfs:
            args.do_not_output_histograms = False
            args.do_not_output_matrix = False

        if args.kr is not None:
            args.inverse_kmer_reference = args.kr
            args.kmer_reference = args.kr

        compare(
            fq_files=args.fq_files, gzipped=args.gzipped,
            no_kmer_filtering=args.no_kmer_filtering,
            coverage_cutoff=args.coverage_cutoff,
            qual_filter=args.qual_score, cpus=args.cpus,
            output_matrix=args.do_not_output_matrix,
            output_histogram=args.do_not_output_histograms,
            output_prefix=args.output_prefix,
            kmer_reference=args.kmer_reference,
            inverse_kmer_reference=args.inverse_kmer_reference,
            pairwise_kmer_filtering=args.pairwise_kmer_filter,
            rapid_mode=args.rapid_mode,
            include_ard_comparison=args.include_ard_comparison,
            jf_input=args.jf_input,
            clean_files=args.keep_files
        )

    elif args.subparser_name == "plot":
        plot_output(input_file=args.input, output_prefix=args.output_prefix)

    elif args.subparser_name == "update_mlst":
        update_mlst_resources()

    else:
        raise AttributeError("{0} is not a valid sub-command".format(args.subparser_name))


if __name__ == "__main__":
    main()
