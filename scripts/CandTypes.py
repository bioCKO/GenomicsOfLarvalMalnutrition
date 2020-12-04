import sys
from collections import defaultdict as d
from optparse import OptionParser, OptionGroup

#Author: Martin Kapun

#########################################################   HELP   #########################################################################
usage="python %prog --glmm cand.glmm --fet cand.fet > output.cand"
parser = OptionParser(usage=usage)
group=OptionGroup(parser,
"""
H E L P :
_________

Annotate candidates according to statistical approach and allele frequency patterns
""")
#########################################################   CODE   #########################################################################

parser.add_option("--glmm", dest="glmm", help="glmm cand list")
parser.add_option("--fet", dest="fet", help="fet cand list")

(options, args) = parser.parse_args()
parser.add_option_group(group)


names=["GLMM","FET"]
files=[options.glmm,options.fet]
CAND=d(lambda: d(lambda: d(list)))
for i in range(len(names)):
    for l in open(files[i],"r"):
        a=l.rstrip().split("\t")
        pops=a[3:]
        CON=[float(x) for x in pops[:6]]
        SEL=[float(x) for x in pops[6:]]
        Test=[1 for x in SEL if x>=0.75]
        Test2=[1 for x in SEL if x>0.9 or x<0.1]
        if sum(Test)>=5:
            CAND[a[0]][int(a[1])]["type"]="High"
            CAND[a[0]][int(a[1])]["Stat"].append(names[i])
        elif Test2 == 6 and sum(SEL)/6 > 0.25 and sum(SEL)<0.75:
            CAND[a[0]][int(a[1])]["type"]="Mid"
            CAND[a[0]][int(a[1])]["Stat"].append(names[i])
        else:
            CAND[a[0]][int(a[1])]["type"]="ambiguous"
            CAND[a[0]][int(a[1])]["Stat"].append(names[i])

for C,V in sorted(CAND.items()):
    for P,V1 in sorted(V.items()):
        print C+"/t"+str(P)+"\t"+V1["type"]+"\t"+",".join(V1["Stat"])
