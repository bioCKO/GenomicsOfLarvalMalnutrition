# Bioinformatics pipeline of Kawecki et al. 2021
The bioinformatics pipeline for the generation and analyses of genomics data in Kawecki et al. 2021.

## A) trim, map and re-align around InDels

We closely followed the bioinformatics pipeline of the data analyses in Kapun *et al.* 2020 [(doi://10.1093/molbev/msaa120)](https://doi.org/10.1093/molbev/msaa120). See the [documentation](https://github.com/capoony/DrosEU_pipeline/blob/master/README.md) for more details

## B) Merging BAM files and joint SNP calling

We closely followed the bioinformatics pipeline of the data analyses in Kapun *et al.* 2020 [(doi://10.1093/molbev/msaa120)](https://doi.org/10.1093/molbev/msaa120). See the [documentation](https://github.com/capoony/DrosEU_pipeline/blob/master/README.md) of the DrosEU bioinformatics pipeline and [documentation](https://github.com/capoony/PoolSNP/blob/master/README.md) of PoolSNP for details. Commands that differ from the standard pipeline and additional scripts are listed below.

```bash
sh /GitHub/PoolSNP/PoolSNP.sh  \
mpileup=Input.mpileup.gz \
reference=Dmel_6.04_hologenome_v2.fasta  \
names=C1,C2,C3,C4,C5,C6,S1,S2,S3,S4,S5,S6 \
max-cov=0.95 \
min-cov=10 \
min-count=5 \
min-freq=0.002 \
miss-frac=0.1 \
jobs=22 \
badsites \
base-quality=15 \
output=SNPdata
```

## C) Calculation of unbiased population genetics estimators Tajima's *pi*, Watterson's *Theta* and Tajima's *D*

We closely followed the bioinformatics pipeline of the data analyses in Kapun *et al.* 2020 [(doi://10.1093/molbev/msaa120)](https://doi.org/10.1093/molbev/msaa120). See the [documentation](https://github.com/capoony/DrosEU_pipeline/blob/master/README.md) for more details

## D) Statistical identification of candidate SNPs based on (i) Fisher's exact test (FET) and (ii) Generalized linear mixed models with a binomial error structure (GLMM)


### 1) calculate GLMM using [GNU parallel](https://www.gnu.org/software/parallel/)

```bash
gunzip -c SNPdata.sync.gz  \
| parallel \
--jobs 22 \
--pipe \
--no-notice \
-k --cat python2.7 /scripts/glmm.py \
--input {} \
--order 1,2,3,4,5,6+7,8,9,10,11,12-1,2,3,7,8,9+4,5,6,10,11,12-1,3,5,7,9,11+2,4,6,8,10,12-3,4,5,9,10,11+1,2,6,7,8,12-2,3,4,8,9,10+1,5,6,7,11,12 \
--group 0,0,0,0,0,0,1,1,1,1,1,1 \
--pop 0,1,2,3,4,5,6,7,8,9,10,11 \
| gzip > SNPdata.glmm.gz
```

### 2) subsample the SNP dataset in SYNC file format to an even 30x read-depth.

See the [documentation](https://github.com/capoony/DrosEU_pipeline/blob/master/README.md) of the DrosEU bioinformatics pipeline for more details

### 3) calculate FET of dataset with normalized read-depth using [GNU parallel](https://www.gnu.org/software/parallel/)

```bash
gunzip -c SNPdata-30x.sync.gz  \
| parallel \
--jobs 22 \
--pipe \
--no-notice \
-k --cat python2.7 /scripts/fet.py \
--input {} \
--order 1,2,3,4,5,6+7,8,9,10,11,12-1,2,3,7,8,9+4,5,6,10,11,12-1,3,5,7,9,11+2,4,6,8,10,12-3,4,5,9,10,11+1,2,6,7,8,12-2,3,4,8,9,10+1,5,6,7,11,12 \
| gzip > SNPdata-30x.fet.gz
```
