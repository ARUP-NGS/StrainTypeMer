__author__ = 'ksimmon'

import jellyfish
import numpy as np

rf = "/home/ksimmon/data/strain_typing/BactieraRun4_22_2014/staph/jellyfish_31_clean/dump_2/" + \
     "recount/merged_db/Staph_merged_q25_31bp.fa"


jf_base = "/home/ksimmon/data/strain_typing/BactieraRun4_22_2014/staph/jellyfish_31_clean/"
db_files = ["S1_S12_clean_q25_31bp.jf",
            "S2_S13_clean_q25_31bp.jf",
            "S3_S14_clean_q25_31bp.jf",
            "S4_S15_clean_q25_31bp.jf",
            "S5_S16_clean_q25_31bp.jf",
            "S6_S17_clean_q25_31bp.jf",
            "S7_S18_clean_q25_31bp.jf",
            "SC_S19_clean_q25_31bp.jf"]

rf = "/home/ksimmon/data/strain_typing/BactieraRun4_22_2014/enterococcus/jellyfish_31_clean/dump_2/recount/merged/entro_merged_q25.fa"

jf_base = "/home/ksimmon/data/strain_typing/BactieraRun4_22_2014/enterococcus/jellyfish_31_clean/"
db_files = ["E1_S20_clean_q25_31bp.jf",
            "E2_S21_clean_q25_31bp.jf",
            "E3_S22_clean_q25_31bp.jf",
            "E4_S23_clean_q25_31bp.jf",
            "EC_S24_clean_q25_31bp.jf"]

query_dbs = []
for f in db_files:
    query_dbs.append(jellyfish.QueryMerFile(jf_base + f))

count_good = 0
count_bad = 0
count_good_singletons = [0,0,0,0,0,0,0,0]
count_bad_singletons  = [0,0,0,0,0,0,0,0]
wf_good = open(rf[:-3] + "_remove_low_counts.txt", "w")
wf_bad = open(rf[:-3] + "_low_counts.txt", "w")
for line in open(rf):
    db_code = line.split("\t")[1].strip()
    kmer = line.split("\t")[0]
    mer = jellyfish.MerDNA(kmer)
    mer.canonicalize()
    answers = []
    for query_db in query_dbs:

        answers.append(str(query_db[mer]))

    v = np.array( answers)
    v[v ==''] = 0.0
    v = v.astype(np.int)

    np.count_nonzero(v)

    if v.max() > 5:
        if np.count_nonzero(v) == 1:
            count_good_singletons[np.nonzero(v)[0][0]] += 1
            #print count_good_singletons
        #wf_good.write(str(mer) + "\t" + db_code + "\t" +"\t".join(answers) + "\n")
        wf_good.write("{0}\t{1}\t{2}\n".format(mer,db_code,"\t".join(answers)))
        count_good += 1
    else:
        if np.count_nonzero(v) == 1:
            count_bad_singletons[np.nonzero(v)[0][0]] += 1
            #print count_bad_singletons
        #wf_bad.write(str(mer) + "\t" + db_code  + "\t" +"\t".join(answers) + "\n")
        wf_bad.write("{0}\t{1}\t{2}\n".format(mer,db_code,"\t".join(answers)))
        count_bad += 1

print "{0}\tgood\n{1}\tbad\n{2}\ttotal".format(count_good,count_bad,count_good+count_bad)
print count_good_singletons
print count_bad_singletons