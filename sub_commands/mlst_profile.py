import jellyfish
import os
from collections import namedtuple, OrderedDict
from Bio import SeqIO
import itertools
import sys


kmer_size = 31

def load_ST(species):
    ST_list = {}
    matching_sequences = {
        "Oxf_gltA"  : [],
        "Oxf_gyrB"  : [],
        "Oxf_gdhB"  : [],
        "Oxf_recA"  : [],
        "Oxf_cpn60" : [],
        "Oxf_gpi"   : [],
        "Oxf_rpoD"  : [],
    }

    _dict = {}

    if species is 'abaumannii':
        with open("/Users/ksimmon/Git/StrainTypeMer/resources/abaumannii/abaumannii_Oxf_ST.txt") as f:
            lines = f.readlines()
            header = lines[0].strip().split("\t")
            print header
            ST = namedtuple('ST', header)
            for line in lines[1:]:
                line = line.strip("\n").split("\t")
                ST_list.update({line[0] : ST(*line)})
                #print ":".join(line[1:-2])
                _dict.update({ ":".join(line[1:-2]) : line[0] })
    return ST_list, _dict


def screen_oxford_sequences(qf):
    base = "/Users/ksimmon/Git/StrainTypeMer/resources/abaumannii/"
    fastas = ["Oxf_gltA.tfa", "Oxf_gyrB.tfa", "Oxf_gdhB.tfa", "Oxf_recA.tfa", "Oxf_cpn60.tfa", "Oxf_gpi.tfa",
              "Oxf_rpoD.tfa"]
    matching_sequences = OrderedDict()
    matching_sequences.update({"Oxf_gltA"  : []})
    matching_sequences.update({"Oxf_gyrB"  : []})
    matching_sequences.update({"Oxf_gdhB"  : []})
    matching_sequences.update({"Oxf_recA"  : []})
    matching_sequences.update({"Oxf_cpn60" : []})
    matching_sequences.update({"Oxf_gpi"   : []})
    matching_sequences.update({"Oxf_rpoD"  : []})
    for fasta in fastas:
        for i, s in enumerate(SeqIO.parse(base + fasta, "fasta")):
            for j in range(0, len(s) - kmer_size + 1):
                kmer = s.seq[j : j+31]
                mer = jellyfish.MerDNA(str(kmer))
                mer.canonicalize()
                if qf[mer] == 0: # break loop if not exact match
                    break
            else:
                gene = "_".join(s.name.split("_")[0:-1])
                seq_num = s.name.split("_")[-1]
                matching_sequences[gene].append(seq_num)

    return [":".join(t) for t in list(itertools.product(*matching_sequences.values()))]


def mlst_profile(jf_file, species):
    _l, st_dict = load_ST(species)


    qf = jellyfish.QueryMerFile(jf_file)

    if species is 'abaumannii':
        matching_sequences = screen_oxford_sequences(qf)

        for profile in matching_sequences:
            if profile in st_dict:
                mat = _l[st_dict[profile]],
                print ",".join(list(mat[0])[:-2])

            else:
                print "NO match: {0}".format(profile)









def main():
    base_jf = "/Users/ksimmon/Box Sync/ARUP/2x150/"

    jf_files = ["A1A_Acinetobacter_baumannii_SRR1314618_all_jellyfish_31.jf",
                "A1B_Acinetobacter_baumannii_SRR1314619_all_jellyfish_31.jf",
                "A2A_Acinetobacter_baumannii_SRR1314620_all_jellyfish_31.jf",
                "A2B_Acinetobacter_baumannii_SRR1314621_all_jellyfish_31.jf",
                "A3A_Acinetobacter_baumannii_SRR1314622_all_jellyfish_31.jf",
                "A3B_Acinetobacter_baumannii_SRR1314623_all_jellyfish_31.jf",
                "A5A_Acinetobacter_baumannii_SRR1314624_all_jellyfish_31.jf",
                "A5B_Acinetobacter_baumannii_SRR1314625_all_jellyfish_31.jf",
                "A6A_Acinetobacter_baumannii_SRR1314626_all_jellyfish_31.jf",
                "A6B_Acinetobacter_baumannii_SRR1314627_all_jellyfish_31.jf",
                "A7A_Acinetobacter_baumannii_SRR1314628_all_jellyfish_31.jf",
                "A7B_Acinetobacter_baumannii_SRR1314629_all_jellyfish_31.jf",
                "A8A_Acinetobacter_baumannii_SRR1314630_all_jellyfish_31.jf",
                "A8B_Acinetobacter_baumannii_SRR1314631_all_jellyfish_31.jf",
                "A9A_Acinetobacter_baumannii_SRR1314632_all_jellyfish_31.jf",
                "A9B_Acinetobacter_baumannii_SRR1314633_all_jellyfish_31.jf",
                "A10A_Acinetobacter_baumannii_SRR1314604_all_jellyfish_31.jf",
                "A10B_Acinetobacter_baumannii_SRR1314605_all_jellyfish_31.jf",
                "A11A_Acinetobacter_baumannii_SRR1314606_all_jellyfish_31.jf",
                "A11B_Acinetobacter_baumannii_SRR1314607_all_jellyfish_31.jf",
                "A12A_Acinetobacter_baumannii_SRR1314608_all_jellyfish_31.jf",
                "A12B_Acinetobacter_baumannii_SRR1314609_all_jellyfish_31.jf",
                "A13A_Acinetobacter_baumannii_SRR1314610_all_jellyfish_31.jf",
                "A13B_Acinetobacter_baumannii_SRR1314611_all_jellyfish_31.jf",
                "A14A_Acinetobacter_baumannii_SRR1314612_all_jellyfish_31.jf",
                "A14B_Acinetobacter_baumannii_SRR1314613_all_jellyfish_31.jf",
                "A15A_Acinetobacter_baumannii_SRR1314614_all_jellyfish_31.jf",
                "A15B_Acinetobacter_baumannii_SRR1314615_all_jellyfish_31.jf",
                "A16A_Acinetobacter_baumannii_SRR1314616_all_jellyfish_31.jf",
                "A16B_Acinetobacter_baumannii_SRR1314617_all_jellyfish_31.jf",]


    species = "abaumannii"
    for f in jf_files:
        sys.stdout.write("{0},".format(f.split("_")[0]))
        mlst_profile(base_jf + f, species)


if __name__ == "__main__":
    main()
