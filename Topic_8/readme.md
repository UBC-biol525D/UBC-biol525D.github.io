---
title: "Topics 8: Population genomics - Structure"
permalink: /Topic_8/
topickey: 8
topictitle: "Population genomics"
---

## Accompanying material
* [Slides](./Topic 8mod.pdf)

## Daily Questions:
1. Whats the relasionship between QUAL and GQ? 
2. What effect does missing data have on a PCA? What about linkage?
3. What distinguishes PCA from Structure-based analyses?

### NOTE:
* If you didn't complete creating _Chinook_GWAS.vcf.gz_ in Topic 7, you can copy it to ~/vcf from /mnt/data/vcf


Last topic we called variants across the three chromosomes. If you look at the VCF, you'll notice there are a lot of sites only genotyped in a small subset of the samples. This can happen with lower overall read depth (in this case this is whole genome sequencing at ~7X depth), but can be due to other factors like divergence been sample and reference. We also have indels, and SNPs with more than two alleles. Many programs strictly require biallelic sites so lets first filter the VCF to a smaller set of usable sites.
We're going to use _bcftools_ a very fast program for processing and filtering VCF files. Here are what we want to filter for:
* At least 80/100 samples genotyped (which is 160 alleles since they're diploid).
* Variant call quality (QUAL) > 30 - a 99.9% chance there is a variant at the site given genotype calls across samples
* Only one alternate allele.
* No indels.
* At least 2 copies of the alternate allele

Lets do a quick check to see how many SNPs we have - so we can keep (loose) track of how many SNPs are removed in each filtering step


```bash
#for BCFtools to work as quickly as it does, VCFs need to be converted to binary format and indexed (for quick referencing)
bgzip Chinook_GWAS.chr_1.vcf
tabix Chinook_GWAS.chr_1.vcf.gz #be careful, gzip also produces .gz suffixes, but won't work with bcftools!

bcftools  view \
	-c 2 \
	-i 'INFO/AN >= 160 && QUAL > 30' \
	-m 2 \
	-M 2 \
	-v snps \
	Chinook_GWAS.chr_1.vcf.gz \
	-O z > Chinook_GWAS.filtered.vcf.gz

#Again, lets index the vcf file for future use
tabix -p vcf Chinook_GWAS.filtered.vcf.gz
```

##Coding challenge
* How many sites remain in the filtered VCF? How many were removed in each step? How many in the chromosome chr_1? Don't forget about grep/zgrep!

A common first pass analysis is to use structure to look at clustering in your data. Admixture is similar to STRUCTURE but orders of magnitude faster. We're going use that, but before that we have to convert our VCF to the bed format. We're going to use plink to do that. Plink is a large set of tools for manipulating genetic data, running GWAS and calculating various stats. It's geared towards human data, so sometimes you have to force it to work with non-human data. For example, it assumes you have human chromosomes (eg 23 of them) and will complain if it doesn't see them.


```bash
cd ~/
mkdir analysis

#we had a bug in our pipeline that appended some extra characters to the beginning of sample names - lets fix this first
#we could use a specialty software like bedtools reheader, but lets just use basic bash commands

zcat vcf/Chinook_GWAS.filtered.vcf.gz | sed 's/-e Chinook/Chinook/g' | bgzip > vcf/Chinook_GWAS.filtered.fixedsamps.vcf.gz

/mnt/bin/plink --make-bed \
	--vcf vcf/Chinook_GWAS.filtered.fixedsamps.vcf.gz \
	--out vcf/Chinook_GWAS_fiiltered_fixedsamps \
	--set-missing-var-ids @:# \
	--double-id \
	--allow-extra-chr
```
This produces files with the suffix .nosex, .log, .fam, .bim, .bed. We can use these in Admixture.

NOTE: When using admixture you should filter your VCF for linkage (i.e. remove highly linked sites). We're going to do this later during the PCA step, so for now we're using our whole set. If you can't filter for linkage, subsetting the site also helps (i.e. selecting every 10th site).

```bash 
/mnt/bin/admixture vcf/Chinook_GWAS_fiiltered_fixedsamps.bed 2
```
Uh oh that doesn't work, it produces this error message.
```bash
****                   ADMIXTURE Version 1.3.0                  ****
****                    Copyright 2008-2015                     ****
****           David Alexander, Suyash Shringarpure,            ****
****                John  Novembre, Ken Lange                   ****
****                                                            ****
****                 Please cite our paper!                     ****
****   Information at www.genetics.ucla.edu/software/admixture  ****

Random seed: 43
Point estimation method: Block relaxation algorithm
Convergence acceleration algorithm: QuasiNewton, 3 secant conditions
Point estimation will terminate when objective function delta < 0.0001
Estimation of standard errors disabled; will compute point estimates only.
Invalid chromosome code!  Use integers.
```
Our chromosomes are named chr_1 and chr_2, not integers like admixture is expecting. Like many programs, this is coded for human data where chromosomes are known and numbered clearly. We need to rename the chromosome column of the vcf so that they're integers. In this case, that means removing "HanXRQChr" from any line that starts with that, although it would depend on how your chromosomes are named.

```bash
zcat vcf/Chinook_GWAS.filtered.fixedsamps.vcf.gz |\
	sed s/^chr_//g |\
	gzip > vcf/Chinook_GWAS.filtered.fixedsamps.numericChr.vcf.gz
	
/mnt/bin/plink --make-bed \
	--vcf vcf/Chinook_GWAS.filtered.fixedsamps.numericChr.vcf.gz \
	--out vcf/Chinook_GWAS.filtered.fixedsamps.numericChr \
	--set-missing-var-ids @:# \
	--double-id \
	--allow-extra-chr

/mnt/bin/admixture --cv vcf/Chinook_GWAS.filtered.fixedsamps.numericChr.bed 2
```
This works! With 100 samples and ~31000 SNPs it finishes relatively quickly. We only ran it for one value of K (2) but we should also test different K values and select the best K value. Common practice is to run from 1 to number of populations (10). Lets do that but skip some in between.
```bash 
for K in 1 2 3 10; \
do /mnt/bin/admixture --cv vcf/Chinook_GWAS.filtered.fixedsamps.numericChr.bed $K |\
tee vcf/Chinook_GWAS.filtered.fixedsamps.numericChr.${K}.out; \
done
#NOTE: "tee" takes the output of a command and saves it to a file, while 
# also printing letting it print to the screen. So we can watch the progress while also 
# saving that output. 

#Now move all the output files to the analysis directory
mv Chinook_GWAS.filtered.fixedsamps.numericChr* analysis/
```
The best K value for Admixture is typically the K value with the lowest cross-validation (CV) error. The CV error are in the .out files we saved. One easy way to look at all those scores is to print all the .out files and then keep only the lines that include "CV" using grep. 

```
cat analysis/*out | grep CV
CV error (K=10): 0.66060
CV error (K=1): 0.40978
CV error (K=2): 0.42452
CV error (K=3): 0.45008
```
This shows that the lowest CV error is with K=1, with K=2 as our second choice. To see how this lines up lets look at the .Q file, which shows group assignment for each sample. The Q file doesn't include sample names so we can put those together using "paste. One way to asses if there is meaningful population structure is to check whether each population grouping has individuals of major ancestry.

```bash
paste <(cut -d" " -f1 vcf/Chinook_GWAS_fiiltered_fixedsamps.fam) analysis/Chinook_GWAS.filtered.fixedsamps.numericChr.2.Q
Chinook.p1.i0	0.148160 0.851840
Chinook.p1.i1	0.000010 0.999990
Chinook.p1.i2	0.000010 0.999990
Chinook.p1.i3	0.000010 0.999990
Chinook.p1.i4	0.000010 0.999990
Chinook.p1.i5	0.000010 0.999990
Chinook.p1.i6	0.000010 0.999990
Chinook.p1.i7	0.000010 0.999990
Chinook.p1.i8	0.000010 0.999990
Chinook.p1.i9	0.000010 0.999990
Chinook.p10.i0	0.999990 0.000010
Chinook.p10.i1	0.999990 0.000010
Chinook.p10.i2	0.999990 0.000010
Chinook.p10.i3	0.999990 0.000010
Chinook.p10.i4	0.999990 0.000010
Chinook.p10.i5	0.999990 0.000010
Chinook.p10.i6	0.821645 0.178355
Chinook.p10.i7	0.999990 0.000010
Chinook.p10.i8	0.999990 0.000010
Chinook.p10.i9	0.999990 0.000010

```
While K=1 has the lowest cross validation error, clearly we can see that individuals from population 1 and population 10 (furthest apart in sampling space) beling to different groupings. 

We're going to plot these results, but before we leave the command line, lets also run another couple of analyses. Lets run a *PCA*. This is a very nice model free approach to visualizing population structure, and is a nice complement to model based structure analyses, like admixture.

We can do this with just one line of code in plink.

```bash
plink --bfile vcf/Chinook_GWAS_fiiltered_fixedsamps --pca --allow-extra-chr --out analysis/Chinook_GWAS_fiiltered_fixedsamps
```

Sweet. Before we get to plotting in R, let's see if we can figure out any regions of the genome that might be under selection in our Chinook salmon. We're going to take two approaches - an Fst scan and Genome-wide association of phenotypes with individual genotypes at each locus (GWA). 



