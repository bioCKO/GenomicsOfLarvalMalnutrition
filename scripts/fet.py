import sys
from collections import defaultdict as d
from collections import Counter
from scipy import stats
from math import asin
from optparse import OptionParser, OptionGroup

#Author: Martin Kapun

#########################################################   HELP   #########################################################################
usage="python %prog --input input.sync --order 1,2,3,4,5,6+7,8,9,10,11,12-1,2,3,7,8,9+4,5,6,10,11,12 > output.fet"
parser = OptionParser(usage=usage)
group=OptionGroup(parser,
"""
H E L P :
_________

Calculate Fisher Exact Tests based on averaged allele counts across all populations within a group, where populations are separated by a ',' and groups by a '+'. Multiple (permuted) combinations of populations across the two groups can be tested. These combinations are separated by a "-". The script will calculate p-values for each of these combinations in the correpsonding order.
""")
#########################################################   CODE   #########################################################################

parser.add_option("--input", dest="sync", help="sync file with all SNPs")
parser.add_option("--norm", dest="norm", help="normalize a la Tad",default="No")
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

def af(x,m):
    y=nuc(x)
    return y.count(m)/float(len(y))

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

def median(x):
    mid =int(len(x)/2)
    sort=sorted(x)
    if len(x)==0:
        return None
    if len(x)%2==0:
        lower=sort[mid-1]
        upper=sort[mid]
        return (float(lower)+float(upper))/2.0
    else:
        return sort[mid]

def fet(a,A,z):
    '''summarize the allele count for alleles a and A for the pairing according to z and perform FST'''
    mac=sum([a[i-1] for i in map(int,z.split("+")[0].split(","))])
    mas=sum([a[i-1] for i in map(int,z.split("+")[1].split(","))])
    mic=sum([A[i-1] for i in map(int,z.split("+")[0].split(","))])
    mis=sum([A[i-1] for i in map(int,z.split("+")[1].split(","))])
    oddsratio, pvalue = stats.fisher_exact([[mac, mas], [mic, mis]],alternative='two-sided')
    return pvalue

def normalize(a,A):
    ''' normalization a la Tad'''
    t=[A[i]+a[i] for i in range(len(a))]
    tmed=median(t)

    tlist=[]
    for i in range(len(a)):
        if t[i]==0:
            tlist.append(0)
        else:
            tlist.append(tmed/t[i])

    an=[round(a[i]*tlist[i]) for i in range(len(a)) ]
    An=[round(A[i]*tlist[i]) for i in range(len(a)) ]
    return an,An


pairings=options.order.split("-")

data=open(options.sync,"r")


for l in data:
    a=l.split()
    if len(a)<3:
        continue
    c,p,r=a[:3]
    lines=a[3:]

    ## remove SNPs with missing data
    if "0:0:0:0:0:0" in lines:
        continue

    alleles=major_alleles(lines)

    if alleles=="na":
        continue
    ## counts for alleles 0 and 1
    a=[nuc(x).count(alleles[0]) for x in lines]
    A=[nuc(x).count(alleles[1]) for x in lines]

    if options.norm!="No":
        # normalize the counts
        an,An=normalize(a,A)

    #print a,an,A,An

    Plist=[]
    for j in range(len(pairings)):
        #print pairings[j]
        #Plist.append(fet(an,An,pairings[j]))
        if options.norm!="No":
            Plist.append(fet(an,An,pairings[j]))
        else:
            Plist.append(fet(a,A,pairings[j]))
    Plist.append(min(Plist[1:]))
    #print "\t".join([c,p,r])+"\t"+"\t".join(map(str,Plist))
    print l.rstrip()+"\t"+"\t".join(map(str,Plist))
