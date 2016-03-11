#!/home/ksimmon/anaconda/bin/python


"""
This script will take as input jf databases from strains and perform a comparison
The input of strains should be provided should be provided as a comma separated list.
assume sample name is before underscore
"""
import sys
from sub_commands.compare_strains import compare_strains
from sub_commands.update_mlst_resources import update_mlst_resources
from sub_commands.compare_fastqs import compare_fastqs
import argparse
import matplotlib
matplotlib.use("Agg")
#jellyfish_path = "/usr/local/bin/jellyfish"


########################################################################################################################
########################################################################################################################

def arguments():
    """
    Parses user supplied arguments for StrainTypeMer subcommands
    :return: named tuple of arguments
    """
    parser = argparse.ArgumentParser(prog="StrainTypeMer", formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description="the script take multiple jellyfish counts and compares the strains")

    subparsers = parser.add_subparsers(help='sub-command help', dest='subparser_name', title="subcommands")

    parser_query = subparsers.add_parser('compare', help='query an NM or GENE in the transcript db', )

    parser_query.add_argument("-c", "--cutoff",
                              help="The number of kmers to analyze [Default = None (e.g analyze all kmers)]", type=int,
                              default=None)

    parser_query.add_argument("--coverage_cutoff",
                              help="percent of genome coverage to set kmer filters [DEFAULT .20 if coverage is 30 " +
                                   "[(30 * .20) = 6] kmers with a count < 5 will be ignored for cooresponding strain",
                              type=float, default=.15)

    parser_query.add_argument("-t", "--cpus",
                              help="The number of cpus to use when counting kmers in strains [Default == 2]", type=int,
                              default=2)

    parser_query.add_argument("--no_kmer_filtering", help="Do not filter kmers based on coverage", action="store_true",
                              default=False)

    parser_query.add_argument("-k", "--kmer_reference",
                              help="supplement use kmer reference set for comparison (e.g. plasmid, core genome, " +
                                   "pan genome kmers)", type=str, default=None, )

    parser_query.add_argument("-r", "--inverse_kmer_reference",
                              help="Use kmers to create similarity matrix not in supplied reference set (e.g. " +
                                   "plasmid, core genome, pan genome kmers)",
                              type=str, default=None, )

    parser_query.add_argument("--do_not_output_histograms",
                              help="This will prevent the output of the PDF files containing the histograms",
                              default=True, action="store_false")

    parser_query.add_argument("--do_not_output_matrix",
                              help="This will prevent the output of the PDF files containing the matrix", default=True,
                              action="store_false")

    parser_query.add_argument("--no_pdfs", help="Output will only goto stdout", default=False, action="store_true")

    parser_query.add_argument("-o", "--output_prefix", help="appends a prefix to the output files", default="")

    parser_query.add_argument('jf_files', nargs='+', help='jellyfish files for each strain')

    subparsers.add_parser('update_mlst_resources', help='update mlst resources', )

    parser_fastq = subparsers.add_parser('compare_fastqs', help='query an NM or GENE in the transcript db', )

    parser_fastq.add_argument("-c", "--cutoff",
                              help="The number of kmers to analyze [Default = None (e.g analyze all kmers)]", type=int,
                              default=None)

    parser_fastq.add_argument("--coverage_cutoff",
                              help="percent of genome coverage to set kmer filters [DEFAULT .20 if coverage is 30 " +
                                   "[(30 * .20) = 6] kmers with a count < 5 will be ignored for cooresponding strain",
                              type=float, default=.15)

    parser_fastq.add_argument("-t", "--cpus",
                              help="The number of cpus to use when counting kmers in strains [Default == 2]", type=int,
                              default=2)

    parser_fastq.add_argument("--no_kmer_filtering", help="Do not filter kmers based on coverage", action="store_true",
                              default=False)

    parser_fastq.add_argument("--qual_score", help="the phred score to filter bases", default=0, type=int)

    parser_fastq.add_argument("-k", "--kmer_reference",
                              help="supplement use kmer reference set for comparison (e.g. plasmid, core genome, " +
                                   "pan genome kmers)", type=str, default=None, )

    parser_fastq.add_argument("-r", "--inverse_kmer_reference",
                              help="Use kmers to create similarity matrix not in supplied reference set (e.g. " +
                                   "plasmid, core genome, pan genome kmers)",
                              type=str, default=None, )

    parser_fastq.add_argument("--do_not_output_histograms",
                              help="This will prevent the output of the PDF files containing the histograms",
                              default=True, action="store_false")

    parser_fastq.add_argument("--do_not_output_matrix",
                              help="This will prevent the output of the PDF files containing the matrix", default=True,
                              action="store_false")

    parser_fastq.add_argument("--no_pdfs", help="Output will only go to stdout", default=False, action="store_true")

    parser_fastq.add_argument("-o", "--output_prefix", help="appends a prefix to the output files", default="")

    parser_fastq.add_argument("-gz", "--gzipped", help="flag to indicate fastq_files are gzipped [Default False]",
                              action="store_true", default=False)

    parser_fastq.add_argument('fq_files', nargs='+',
                              help='fastq files for each strain (fq1 OR fq1;label OR fq1,fq2[label will '
                                                          'be matching string of the two files OR fq1,fq2;label')

    args = parser.parse_args()

    return args
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


def main():
    args = arguments()

    if args.subparser_name == "compare":
        # set output of pdfs to false
        if args.no_pdfs:
            args.do_not_output_histogram = False
            args.do_not_output_matrix = False

        compare_strains(jf_files=args.jf_files,
                        no_kmer_filtering=args.no_kmer_filtering,
                        cutoff=args.cutoff,
                        cpus=args.cpus, coverage_cutoff=args.coverage_cutoff,
                        kmer_reference=args.kmer_reference,
                        inverse_kmer_reference=args.inverse_kmer_reference,
                        output_histogram=args.do_not_output_histograms,
                        output_matrix=args.do_not_output_matrix,
                        output_prefix=args.output_prefix,
                        )

    if args.subparser_name == "compare_fastqs":
        compare_fastqs(fq_files=args.fq_files, gzipped=args.gzipped, coverage_cutoff=args.coverage_cutoff,
                       qual_filter=args.qual_score, cpus=args.cpus)

    elif args.subparser_name == "update_mlst_resources":
        update_mlst_resources()

    else:
        pass


if __name__ == "__main__":

    main()

