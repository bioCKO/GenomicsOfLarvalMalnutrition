import sys
from collections import defaultdict as d
from optparse import OptionParser, OptionGroup

#Author: Martin Kapun

#########################################################   HELP   #########################################################################
usage="python %prog --input input.sync > output.af"
parser = OptionParser(usage=usage)
group=OptionGroup(parser,
"""
H E L P :
_________

Calculate allele frequencies of the alleles more common in the selected population than in the control populations. The script assumes that there is an equal number of selected and control in the input sync file and that the first half is comprised by the control and the second by the selected populations. The output contains the allelic states of the "Selected" and "Control" alleles (separated by a dash) in the third column
""")
#########################################################   CODE   #########################################################################

parser.add_option("--input", dest="input", help="input SYNC file")

(options, args) = parser.parse_args()
parser.add_option_group(group)

def sync2freqh(x):
    ''' convert string in SYNC format to dictionary of freqencies where x is a string in sync format'''
    from collections import defaultdict as d

    nuc=["A","T","C","G"]
    counts=map(int,x.split(":")[:4])
    if sum(counts)==0:
        return ({"A":0.0,"T":0.0,"C":0.0,"G":0.0},0)
    CO=dict(zip(*[nuc,counts]))
    h=d(float)
    for k,v in CO.items():
        h[k]=v/float(sum(CO.values()))

    return h,sum(CO.values())

def all_alleles(v):
    ''' returns most common strings'''
    from collections import Counter
    nuc=v.replace("N","")
    return zip(*Counter(nuc).most_common())[0]

def sync2string(x):
    ''' convert sync format to string of nucleotides  where x is a string in sync format '''
    string=""
    alleles=["A","T","C","G"]
    ah=dict(zip(alleles,map(int,x.split(":")[:4])))
    for k,v in ah.items():
        string+=v*k
    return string

def load_data(x):
    ''' import data either from a gzipped or or uncrompessed file or from STDIN'''
    import gzip
    if x=="-":
        y=sys.stdin
    elif x.endswith(".gz"):
        y=gzip.open(x,"r")
    else:
        y=open(x,"r")
    return y

def find_sel_allele(Con,Sel,M,m):
    Sf,Cf=[],[]
    for i in Sel:
        Sf.append(sync2freqh(i)[0][M])
    for i in Con:
        Cf.append(sync2freqh(i)[0][M])

    if sum(Cf)/float(len(Con)) < sum(Sf)/float(len(Sf)):
        return M,m
    else:
        return m,M

for l in load_data(options.IN):
    a=l.split()
    ID=a[0]+":"+a[1]
    AlleleStrings=""
    pops=a[3:]
    for pop in pops:
        AlleleStrings+=sync2string(pop)
    Ma,ma=all_alleles(AlleleStrings)[:2]
    Ma,ma=find_sel_allele(pops[:len(pops)/2],pops[len(pops)/2:],Ma,ma)
    fl=[]
    for pop in a[3:]:
        FH=sync2freqh(pop)[0]
        if "na" in FH:
            fl.append("NA")
            continue
        fl.append(FH[Ma])
    #print l
    print "\t".join(a[:2])+"\t"+Ma+"/"+ma+"\t"+"\t".join(map(str,fl))
