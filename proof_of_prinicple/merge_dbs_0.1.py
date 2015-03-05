__author__ = 'ksimmon'

import math

import shutil as su
import subprocess

files=["A1_S1_31bp_2_count.jf",
"A2_S2_31bp_2_count.jf",
"A3_S3_31bp_2_count.jf",
"A42_S10_31bp_3_count.jf",
"A4_S3_31bp_1_count.jf",
"A52_S11_31bp_2_count.jf",
"A5_S5_31bp_2_count.jf",
"A6_S6_31bp_2_count.jf",
"A7_S7_31bp_2_count.jf",
"A8_S8_31bp_2_count.jf",
"A9_S9_31bp_2_count.jf",
"AC1_S1_31bp_3_count.jf",
"AC2_S2_31bp_3_count.jf",
"AC3_S3_31bp_3_count.jf",
"AC4_S4_31bp_3_count.jf",
"AC5_S1_31bp_3_count.jf",
"AC6_S2_31bp_3_count.jf",
"AC7_S3_31bp_3_count.jf",
"AC8_S4_31bp_3_count.jf",
"AC9_2_S5_31bp_3_count.jf",
"AC9_S5_31bp_3_count.jf"
]
#kl db mod
kl = 31

_t = 0
for i in range(len(files)):
    offset = int(math.pow(2,i))
    print offset
    _t += offset
    print files[i]
    file_mv = files[i][:-3] + "_{0}.jf".format(offset)
    #subprocess.call("python /home/ksimmon/bin/Taxonomer/utilities/mod_jf2.0_count.py " \
            # + str(kl) + " " + files[i] + " " + str(offset), shell=True)

    #su.move(file[i], file_mv)

print _t