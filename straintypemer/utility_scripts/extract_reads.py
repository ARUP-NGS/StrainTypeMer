

import argparse
import sys
from Bio import SeqIO

parser = argparse.ArgumentParser(description="extract reads from genome reads from kraken results")

parser.add_argument("--kraken_input", help="kraken output", type=argparse.FileType('r'))
parser.add_argument("--fasta_file", help="fasta file to extract sequencegres", type=str)
parser.add_argument('--taxonomy_file',  help='jellyfish files for each strain', type=argparse.FileType('r'))
parser.add_argument('--parent_taxid',  help='parent taxid to extract', type=int)

args = parser.parse_args()
kraken_input = args.kraken_input
fasta_file = args.fasta_file
taxonomy_file = args.taxonomy_file
parent_taxid = args.parent_taxid

# base = "/home/ksimmon/data/strain_typing/MiSeq_run_July_2015/renamed/148-6_3-Acinetobacter-gp3_processed_jellyfish_31/"
# kraken_input = base + "148-6_3-Acinetobacter-gp3_kraken.krk"
# fasta_file = base + "148-6_3-Acinetobacter-gp3_preprocessed_all.fasta"
# taxonomy_file = "/home/ksimmon/reference/strian_typing_resources/kraken_acineto_db/taxonomy/nodes.dmp"
# parent_taxid = 2


tax_tree = {}

for line in taxonomy_file:
    line = line.strip().split("\t")
    child, parent = int(line[0]), int(line[2])
    tax_tree.update({child:parent})

reads_to_extract = set([])
for line in kraken_input:
    line = line.strip().split("\t")
    is_classified, read_name, taxid = line[0] == "C", line[1], int(line[2])
    if is_classified:
        while taxid != 1:
            if taxid == parent_taxid:
                reads_to_extract.add(read_name)
                break

            taxid = tax_tree[taxid]

for s in SeqIO.parse(fasta_file, "fasta"):
    if s.name in reads_to_extract:
        sys.stdout.write(">{0}\n{1}\n".format(s.name,s.seq,))
