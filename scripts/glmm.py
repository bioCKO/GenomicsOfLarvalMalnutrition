
import sys
from collections import defaultdict as d
from collections import Counter
from scipy import stats
from math import asin
from rpy2.robjects import r
import rpy2.robjects as robjects
from optparse import OptionParser, OptionGroup

#Author: Martin Kapun

#########################################################   HELP   #########################################################################
usage="python %prog --input input.sync --order 1,2,3,4,5,6+7,8,9,10,11,12-1,2,3,7,8,9+4,5,6,10,11,12 --group 0,0,0,0,0,0,1,1,1,1,1,1 --pop 0,1,2,3,4,5,6,7,8,9,10,11 > output.glmm"
parser = OptionParser(usage=usage)
group=OptionGroup(parser,
"""
H E L P :
_________

requires the rpy2 package

Calculate Generalized linear mixed models with a binomial error structure based on allele counts across all populations within a group, where populations are separated by a ',' and groups by a '+'. The levels of the group effects need to be provide with the parameter (--group) and ID's of the random factor "population" with (--pop).

Multiple (permuted) combinations of populations across the two groups can be tested. These combinations are separated by a "-". The script will calculate p-values for each of these combinations in the correpsonding order.
""")
#########################################################   CODE   #########################################################################

parser.add_option("--input", dest="sync", help="sync file with all SNPs")
parser.add_option("--pop", dest="pop", help="random effect: population ID")
parser.add_option("--group", dest="group", help="fixed effect: Treatment group")
parser.add_option("--order", dest="order", help="the order of the populations")

(options, args) = parser.parse_args()
parser.add_option_group(group)

def nuc(z):
    al=["A","T","C","G"]
    nu=""
    cov=map(int,z.split(":")[:4])
    for i in range(4):
        nu+=al[i]*int(cov[i])
    return nu

def ac(x,m):
    y=nuc(x)
    return (y.count(m[0]),(y.count(m[1])))


def major_alleles(v):
    allnuc=""
    for po in v:
        allnuc+=nuc(po)
    if allnuc=="":
        return "na"
    toptwo=Counter(allnuc).most_common()[:2]
    if len(toptwo)<2:
        return "na"
    return zip(*toptwo)[0]

def glm(x,y):

    z=[x[int(i)-1] for i in y.split(",")]
    r.assign("A1",robjects.vectors.FloatVector(zip(*z)[0]))
    r.assign("A2",robjects.vectors.FloatVector(zip(*z)[1]))

    r.assign("eff",robjects.vectors.StrVector(effect))
    r.assign("pop",robjects.vectors.IntVector(pop))
    r('gm=suppressMessages(tryCatch(glmer(cbind(A1,A2)~eff+(1|pop),family="binomial"),error=function(e) e))')
    if "error" in r('class(gm)'):
        return "na"
    r('gm2=suppressMessages(tryCatch(glmer(cbind(A1,A2)~(1|pop),family="binomial"),error=function(e) e))')
    if "error" in r('class(gm2)'):
        return "na"
    return list(r('suppressMessages(lrtest(gm,gm2))$"Pr(>Chisq)"[2]'))[0]

pairings=options.order.split("-")
effect=options.group.split(",")
pop=map(int,options.pop.split(","))

data=open(options.sync,"r")

r('suppressMessages(library(lme4))')
r('suppressMessages(library(aod))')
r('suppressMessages(library("lmtest"))')

for l in data:
    a=l.split()
    if len(a)<3:
        continue
    c,p,ro=a[:3]
    lines=a[3:]

    ## remove SNPs with missing data
    if "0:0:0:0:0:0" in lines:
        continue

    major=major_alleles(lines)
    if major=="na":
        continue
    allele=major
    #print allele
    alleles=[ac(i,allele) for i in lines]
    #print alleles
    Plist=[]
    for j in range(len(pairings)):
        Plist.append(glm(alleles,pairings[j]))
    Plist.append(min(Plist[1:]))
    if "na" in Plist:
        continue
    print "\t".join([c,p,ro])+"\t"+"\t".join(map(str,Plist))
    #print l.rstrip()+"\t"+"\t".join(map(str,Plist))
