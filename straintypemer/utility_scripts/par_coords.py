data_file = open("/data2/wash_u_sra/ENTEROCOCCUS/out_mar_16_vre.txt")
process = False
d = {}
tag_list = set([])
for line in data_file:
    line = line.strip()
    if "STRAIN STATS" in line:
        #print(line)
        process = True
    #print(process)
    if process:
        #print line
        if "Strain:" in line:
            strain = line.strip().split(" ")[-1]

        if "ARD GENE:" in line:
            tag = line[ line.find("tag:") + 5 : line.find("Covered:") -1 ]
            # tag = tag.split("-")[0]
            tag_list.add(tag)
            covered = float(line[ line.find("Covered:") + 9 : line.find("%")  ])
            #print(strain, tag, covered)
            if strain in d:
                d[strain].update({tag : covered})
            else:
                d.update({strain : {"mlst" : 0, tag : covered }})
        if "MLST profile" in line:
            try:
                profile = int(line[ line.find("ST:") + 4 :  line.find("\tprofile:") ])
            except:
                profile = 0
            #print profile

            if strain in d:
                d[strain].update({"mlst" : profile})
            else:
                d.update({strain : {"mlst" : profile} })


for tag in tag_list:
    for strain in d.keys():
        #print d
        if tag not in d[strain]:
            d[strain].update({tag : 0.0})
import sys
#wf = sys.stdout #
wf = open("/data2/wash_u_sra/ENTEROCOCCUS/out_mar_16_vre_results.csv", "w")
cols = ['mlst'] + list(tag_list)

rows = {}
for strain in sorted(d.keys()):
    s = ""
    for c in cols:
        s += "{0},".format(d[strain][c])
    if s not in rows:
        rows.update({ s: [strain] })
    else:
        rows[s].append(strain)
    #wf.write(s[:-1] + "{0}\n".format(strain))
wf.write(",".join(cols) + ",Strain\n")

_dict = {}
for k,v in rows.iteritems():

    _dict.update({":".join(v) : [ i for i in enumerate(k.split(","))]})
    print "{0}{1}\n".format(k, ":".join(v))
    wf.write("{0}{1}\n".format(k, ":".join(v)))
wf.close()

from pandas.tools.plotting import andrews_curves
import pandas
import matplotlib.pyplot as plt
from pandas.tools.plotting import parallel_coordinates
rf = pandas.read_csv("/data2/wash_u_sra/ENTEROCOCCUS/out_mar_16_vre_results.csv", sep=",")
#print rf
#fig = plt.figure(num=1, figsize=13,13, )
fig = plt.figure()
ax = fig.add_axes([.1, .3,.8, .3])
c = ["#1f78b4", "#33a02c","#e31a1c", "#ff7f00", "#6a3d9a", "#a6cee3", "#b2df8a","#fb9a99","#fdbf6f", "#cab2d6"   ]
p = parallel_coordinates(rf, 'Strain', axvlines_kwds={'linewidth':.5,'color':'black'},axvlines=False, lw=5, alpha=.9,
                         colors=c, antialiased=True)

#plt.setp(p,linewidth=2)

#print rf.columns
p.legend(loc='upper right', bbox_to_anchor=(1,2.2))
p.set_xticklabels(rf.columns, rotation=90)

# for k, v in p.__dict__.iteritems():
#     print k,v


# print p.get_lines()
# print p.get_lines()[0].__dict__
# print p.get_lines()[0]._label
# p.get_lines()[0]._linewidth = 53
# print p.get_lines()[0]._linewidth
#p.set_xticklabels(rotation=90)
#andrews_curves(rf, 'Strain')
#plt.set_xticklabels(p,rotation=90)
#fig.tight_layout()

for i, l in enumerate(p.get_lines()):
    p.get_lines()[i]._linewidth = len(p.get_lines()[i]._label.split(":")) * 1.0


plt.show()