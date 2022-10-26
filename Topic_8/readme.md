---
title: "Topics 8: Population genomics - Structure"
permalink: /Topic_8/
topickey: 8
topictitle: "Population genomics - Structure"
---

## Accompanying material
* [Slides](https://github.com/UBC-biol525D/UBC-biol525D.github.io/blob/master/Topic_8/PopGen%202021.pdf)


Last topic we called variants across the two chromosomes. If you look at the VCF, you'll notice there are a lot of sites only genotyped in a small subset of the samples. This can happen with lower overall read depth (in this case this is whole genome sequencing at ~8X depth), but can be due to other factors like divergence been sample and reference. We also have indels, and SNPs with more than two alleles. Many programs strictly require biallelic sites so lets first filter the VCF to a smaller set of usable sites.

We're going to use _bcftools_ a very fast program for processing and filtering VCF files. Here are what we want to filter for:
* At least 80/100 samples genotyped (which is 160 alleles since they're diploid).
* Variant call quality (QUAL) > 30 - a 99.9% chance there is a variant at the site given genotype calls across samples
* Only one alternate allele.
* No indels.
* At least 2 copies of the alternate allele

### NOTE:
* We are going to be doing some data hungry analyses, so will be working with a much larger VCF then we made last tutorial. Copy it to from `/mnt/data/vcf/Chinook_GWAS.vcf.gz` to `~/vcf`


```bash
cd ~/vcf
bcftools view \
    -c 2 \
    -i 'INFO/AN >= 160 && QUAL > 30' \
    -m 2 \
    -M 2 \
    -v snps \
    -O z \
    Chinook_GWAS.vcf.gz > Chinook_GWAS_filtered.vcf.gz

#Lets index the vcf file for future use
tabix Chinook_GWAS_filtered.vcf.gz
```

##Coding challenge
* How many sites remain in the filtered VCF? How many were removed? How many on each chromosome? Don't forget about `grep -v` ad `wc -l`!

A common first pass analysis is to use structure to look at clustering in your data. Admixture is similar to STRUCTURE but orders of magnitude faster. We're going use that, but before that we have to convert our VCF to the specialized format. We can do that with Plink - Plink is a large set of tools for manipulating genetic data, running GWAS and calculating various stats. It's geared towards human data, so sometimes you have to force it to work with non-human data. For example, it assumes you have human chromosomes (eg 23 of them) and will complain if it doesn't see them.


```bash
cd ~/
mkdir analysis

#we had a bug in our pipeline that appended some extra characters to the beginning of sample names
#you can check sample names by greping for "#CHROM", which is the first string of the sample header line
zgrep "#CHROM" vcf/Chinook_GWAS_filtered.vcf.gz
#we could use a specialty software like bedtools reheader to fix this, but lets just use basic bash commands

zcat vcf/Chinook_GWAS_filtered.vcf.gz | sed 's/-e Chinook/Chinook/g' | bgzip > vcf/Chinook_GWAS_filtered_fixedsamps.vcf.gz
#the key command here is the sed 's/find/replace/g' , zcat is uncompression to standard out and bgzip is recompressing 

plink=/mnt/software/plink
$plink --make-bed \
	--vcf vcf/Chinook_GWAS_filtered_fixedsamps.vcf.gz \
	--out vcf/Chinook_GWAS_filtered_fixedsamps \
	--set-missing-var-ids @:# \
	--double-id \
	--allow-extra-chr 
```
This produces files with the suffix .nosex, .log, .fam, .bim, .bed. We can use these in Admixture.

NOTE: When inferring patterns of population structure (i.e. admixture/pca) its good practice to filter your VCF for linkage (i.e. remove highly linked sites). If you can't filter for linkage, subsetting the site also helps (i.e. selecting every 10th site). Not only does this make your inference draw from independent data, but it also keeps run times down!

To prune for LD, we'll ask plink to slide across the genome (10 snps at a time),and in windows of 100snps, calculate LD between each snp, removing those with an LD (r2) > .5

```bash
$plink \
	--bfile vcf/Chinook_GWAS_filtered_fixedsamps \
	--indep-pairwise 100 10 0.5 \
	--out vcf/Chinook_GWAS_filtered_fixedsamps \
	--make-founders \
	--allow-extra-chr 
#allow extra chromosome is another way we force plink to work with non-human data (i.e. allow chromosomes that don't have the typical human chr names)

#this produces two files, .in (to include) and .out (to exclude) 
#now lets actually extracts the snps that remain after LD pruning

$plink \
--bfile vcf/Chinook_GWAS_filtered_fixedsamps \
--extract vcf/Chinook_GWAS_filtered_fixedsamps.prune.in \
--make-bed \
--out vcf/Chinook_GWAS_filtered_fixedsamps_LDpruned \
--allow-extra-chr 

#great! we have our files in the necessary format, now lets run admixture
admixture=/mnt/software/admixture
$admixture vcf/Chinook_GWAS_filtered_fixedsamps_LDpruned.bed 2
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
Our chromosomes are named chr_1 and chr_2, not integers like admixture is expecting. Like many programs, this is coded for human data where chromosomes are known and numbered clearly. We need to rename the chromosome column of the vcf so that they're integers. In this case, that means removing "chr_" from any line that starts with that, although it would depend on how your chromosomes are named. We'll also have to redo the LD pruning.

```bash
#numeric chr names
zcat vcf/Chinook_GWAS_filtered_fixedsamps.vcf.gz |\
	sed s/^chr_//g |\
	gzip > vcf/Chinook_GWAS_filtered_fixedsamps_numericChr.vcf.gz
	
#sed s/^chr_//g <- removes "chr_" from the beginning of lines 
	
#make new bed from vcf
$plink --make-bed \
	--vcf vcf/Chinook_GWAS_filtered_fixedsamps_numericChr.vcf.gz \
	--out vcf/Chinook_GWAS_filtered_fixedsamps_numericChr \
	--set-missing-var-ids @:# \
	--double-id \
	--allow-extra-chr
	
#redo LD analysis
$plink \
	--bfile vcf/Chinook_GWAS_filtered_fixedsamps_numericChr \
	--out vcf/Chinook_GWAS_filtered_fixedsamps_numericChr \
	--allow-extra-chr \
	--make-founders \
	--indep-pairwise 100 10 0.5 
	
#make new bed with snps in low LD
$plink \
	--bfile vcf/Chinook_GWAS_filtered_fixedsamps_numericChr \
	--extract vcf/Chinook_GWAS_filtered_fixedsamps_numericChr.prune.in \
	--make-bed \
	--out vcf/Chinook_GWAS_filtered_fixedsamps_numericChr_LDpruned \
	--allow-extra-chr
	
#run admixture
$admixture --cv vcf/Chinook_GWAS_filtered_fixedsamps_numericChr_LDpruned.bed 2 |\
tee analysis/Chinook_GWAS_filtered_fixedsamps_numericChr_LDpruned.2.out; \

#NOTE: "tee" takes the output of a command and saves it to a file, while 
# also printing letting it print to the screen. So we can watch the progress while also 
# saving that output. 
```
This works! With 100 samples and ~31000 SNPs across two chromosomes it finishes in less than a minute. We only ran it for one value of K (2) but we should also test different K values and select the best K value. Common practice is to run from 1 to number of populations (10). Lets do that but skip some in between. This will take a couple minutes.

```bash 
for K in 1 3 4 10; \
do 
$admixture --cv vcf/Chinook_GWAS_filtered_fixedsamps_numericChr_LDpruned.bed $K |\
tee analysis/Chinook_GWAS_filtered_fixedsamps_numericChr_LDpruned.${K}.out; \
done

#now move the admixture output files to the analysis folder
mv Chinook_GWAS_filtered_fixedsamps_numericChr_LDpruned*.P Chinook_GWAS_filtered_fixedsamps_numericChr_LDpruned*.Q analysis/

```
The best K value for Admixture is typically the K value with the lowest cross-validation (CV) error. The CV error are in the .out files we saved. One easy way to look at all those scores is to print all the .out files and then keep only the lines that include "CV" using grep. 

```
cat analysis/Chinook_GWAS_filtered_fixedsamps_numericChr_LDpruned*out | grep CV
CV error (K=1): 0.49626
CV error (K=10): 0.79357
CV error (K=2): 0.51421
CV error (K=3): 0.54330
CV error (K=4): 0.57121
```
This shows that the lowest CV error is with K=1, but actually K=2 is a close second. To see how this lines up lets look at the .Q file, which shows group assignment for each sample. The Q file doesn't include sample names so we can put those together using "paste. One way to asses if there is meaningful population structure is to check whether each population grouping has individuals of major ancestry.

```bash
paste <(cut -d" " -f1 vcf/Chinook_GWAS_filtered_fixedsamps.fam) analysis/Chinook_GWAS_filtered_fixedsamps_numericChr_LDpruned.2.Q

#what was that <( ) notation? 
#this is another cool bash trick called *subprocessing*. 
#Instead of making an intermediate file of sample names to join with our K values...
#...we can use subprocessing to take the output of the cut command and pass that to paste. 

#Chinook.p1.i0	0.999990 0.000010
#Chinook.p1.i1	0.999990 0.000010
#Chinook.p1.i2	0.999990 0.000010
#Chinook.p1.i3	0.999990 0.000010
#Chinook.p1.i4	0.999990 0.000010
#Chinook.p1.i5	0.999990 0.000010
#Chinook.p1.i6	0.999990 0.000010
#Chinook.p1.i7	0.999990 0.000010
#Chinook.p1.i8	0.999990 0.000010
#Chinook.p1.i9	0.999990 0.000010
#Chinook.p10.i0	0.000010 0.999990
#Chinook.p10.i1	0.000010 0.999990
#Chinook.p10.i2	0.000010 0.999990
```
While K=1 has the lowest cross validation error, clearly we can see that individuals from population 1 and population 10 (furthest apart in sampling space) beling to different groupings. 

We're going to plot these results, but before we leave the command line, lets also one more analysis. Lets run a *PCA*. This is a very nice model free approach to visualizing population structure, and is a good complement to model based structure analyses, like admixture.

We can do this with just one line of code in plink.

```bash
$plink \
	--bfile vcf/Chinook_GWAS_filtered_fixedsamps_numericChr \
	--pca \
	--allow-extra-chr \
	--out analysis/Chinook_GWAS_filtered_fixedsamps_numericChr
```

This was on our full dataset, but we learned above it might be a good idea to prune for linkage. Lets redo the analysis, but now with the pruned set

```bash
#we already have an LD pruned bed, so subbing in that input file
$plink \
	--bfile vcf/Chinook_GWAS_filtered_fixedsamps_numericChr_LDpruned \
	--pca \
	--allow-extra-chr \
	--out analysis/Chinook_GWAS_filtered_fixedsamps_numericChr_LDpruned
````



That's probably all we have time for in this tutorial, but next time lets get right into R to visualize our results from the admixture and PCA analyses.
