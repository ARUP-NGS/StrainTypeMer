
def write_line(array_of_lines):
    print array_of_lines[0]
    if array_of_lines[0][-1] == "+":
        for i, j in enumerate(array_of_lines):
            pass
            #wf.write("{0}\t{1}_{2}_{3}_{4}\t{5}\t{6}\n".format("\t".join(j[0:3]), j[4], j[5], str(1+i), j[0], "0", "+"))
            wf.write("{0}\t{1}_{2}\t{3}\t{4}\n".format("\t".join(j[0:3]), j[4], j[5], "0", "+"))
    else:
        for i, j in enumerate(array_of_lines):
            wf.write("{0}\t{1}_{2}\t{3}\t{4}\n".format("\t".join(j[0:3]),j[4],j[5], "0", "-"))
            #wf.write("{0}\t{1}_{2}_{3}_{4}\t{5}\t{6}\n".format("\t".join(j[0:3]), j[4], j[5],
             #                                                   str(len(array_of_lines) - i),
                                                              #  j[0], "0", "-"))


"""
This script creates reformatted to BED file, which modifies the EXON numbering native to UCSC genome browser input.

1) Starts the exon count from 1 instead of zero
2) Numbers exons from 5' end of the transcript instead of genome postion.

@input bed file to reformat
@output reformated bed file
"""
__author__ = "Keith Simmon"
__date__ = "April 7, 2015"

import argparse


data_file = "/home/ksimmon/Desktop/4-plex_complete_04_09_2015_coding_exons.bed"
wf = open(".".join(data_file.split(".")[:-1]) + "_formatted_CDS.bed", "w")


line_arr = []
start_nm = ""
first_line = True
for line in open(data_file):
    l = line.strip().split("\t")
    nm = "_".join(l[3].split("_")[0:2])
    annotation_type = l[3].split("_")[2].upper()
    if nm == start_nm:
        line_arr.append(l[0:4] + [nm, annotation_type]  + [l[-1]])
    else:
        if not first_line:
            write_line(line_arr)
        first_line = False
        line_arr = []
        start_nm = nm
        line_arr.append(l[0:4] + [nm, annotation_type] + [l[-1]])
else:
    write_line(line_arr)
