#!/home/ksimmon/anaconda/bin/python


"""
This script will take as input jf databases from strains and perform a comparison
The input of strains should be provided should be provided as a comma separated list.
assume sample name is before underscore
"""
import sys
from sub_commands.compare import compare
from sub_commands.update_mlst_resources import update_mlst_resources
# from sub_commands.compare_fastqs import compare_fastqs
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


    subparsers.add_parser('update_mlst', help='update mlst resources', )

    parser_fastq = subparsers.add_parser('compare', help='compare strains', )

    parser_fastq.add_argument("--coverage_cutoff",
                              help="percent of genome coverage to set kmer filters [DEFAULT .20 if coverage is 30 " +
                                   "[(30 * .20) = 6] kmers with a count < 5 will be ignored for cooresponding strain",
                              type=float, default=.15)

    parser_fastq.add_argument("-t", "--cpus",
                              help="The number of cpus to use when counting kmers in strains [Default == 2]", type=int,
                              default=2)

    parser_fastq.add_argument("--no_kmer_filtering", help="Do not filter kmers based on coverage", action="store_true",
                              default=False)

    parser_fastq.add_argument("-q", "--qual_score", help="the phred score to filter bases", default=0, type=int)

    parser_fastq.add_argument("-k", "--kmer_reference",
                              help="supplement use kmer reference set for comparison (e.g. plasmid, core genome, " +
                                   "pan genome kmers) in fasta_format", type=str, default=None, )

    parser_fastq.add_argument("-r", "--inverse_kmer_reference",
                              help="Use kmers to create similarity matrix not in supplied reference set (e.g. " +
                                   "plasmid, core genome, pan genome kmers) in fastq format",
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

    parser_fastq.add_argument("-pwf", "--pairwise_kmer_filter",
                                help="evaluate non-shared kmers in closely related strains",
                                action="store_true", default=False)

    parser_fastq.add_argument('fq_files', nargs='+',
                              help='fastq files for each strain (fq1 OR fq1;label OR fq1,fq2;label will '
                                                          'be matching string of the two files OR fq1,fq2;label'
                                                          " inclued '[NF]' at the end to prevent kmer filter (useful"
                                   "when adding a reference genome")

    args = parser.parse_args()
    return args
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


def main():
    args = arguments()

    if args.subparser_name == "compare":
        if args.no_pdfs:
            args.do_not_output_histogram = False
            args.do_not_output_matrix = False

        compare(fq_files=args.fq_files, gzipped=args.gzipped,
                no_kmer_filtering=args.no_kmer_filtering,
                coverage_cutoff=args.coverage_cutoff,
                qual_filter=args.qual_score, cpus=args.cpus,
                output_matrix=args.do_not_output_matrix,
                output_prefix=args.output_prefix,
                kmer_reference=args.kmer_reference,
                inverse_kmer_reference=args.inverse_kmer_reference,
                pairwise_kmer_filtering=args.pairwise_kmer_filter,
                )

    elif args.subparser_name == "update_mlst":
        update_mlst_resources()

    else:
        pass


if __name__ == "__main__":

    main()


