
import os
import sys
import argparse
from straintypemer.sub_commands.compare import compare

def arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("csv_file",
                        help="the csv file containing the 'accession number', 'strain id' and 'label' to compare",
                        type=argparse.FileType("rU"))
    parser.add_argument("--delimiter", help="The delimiter value", default=",",)



def parse_csv(csv_file, delimeter):
    files = {}
    for line in csv_file:
        accession, name, label = csv_file.strip().split(delimeter)[0:3]
        files.update({"{0}_{1}".format(accession, name): label})
    return files



def find_files(files_to_find):
    files_to_analyze = {}
    found_files = 0
    base_directory = "/data2/playground/test_structure"
    for directory, _, files in os.walk(base_directory, followlinks=True):
        for f in files:
            name = f.split(".")[0]
            if name in f:
                found_files += 1
                files_to_analyze.update({ name :"{0}:{1}".format(files_to_find, files_to_find[name])})
                sys.stderr.write("found {0}\n".format(name))
    if len(set(files_to_analyze).difference(files_to_find)):
        raise ValueError("the following file was not accession [{0}] and strain [{1} were not found".format(name.split("_")[0], name.split("_")[1]))
    return files_to_analyze





def main():
    args = arguments()
    files_to_find = parse_csv(args.csv_file, args.delimiter)

    files = find_files(files_to_find)




    compare(
        fq_files=list(files.values()), gzipped=True,
        no_kmer_filtering=False,
        coverage_cutoff=0.3,
        qual_filter=0, cpus=10,
        output_matrix=False,
        output_histogram=False,
        output_prefix="test",
        kmer_reference=None,
        inverse_kmer_reference=None,
        pairwise_kmer_filtering=True,
        rapid_mode=False,
        include_ard_comparison=False,
        jf_input=False,
        clean_files=False
    )


if __name__ == '__main__':
    main()



