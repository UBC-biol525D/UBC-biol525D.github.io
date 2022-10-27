---
title: "Topics 9: Population genomics - Selection"
permalink: /Topic_GWAS/
topickey: 9
topictitle: "Population genomics - Selection"
---

Last topic we used bcftools to filter snps, and plink & admixture to investigate population structure. We will plot these results to get a better visual sense of this data, but not quite yet! Before we move to R, we are going to run some analyses to infer potential patterns of selection and adaptation in our dataset. 

Remember, we want to test hypothesis about how our populations of Chinook salmon might have adapted to different temperature environments. We could do this a bunch of different ways, but lets try two. The approach we'll take is to see if we can identify loci that show particularly extreme differentiation between environments of varying temperatures. We'll use the Fst statistic to find these loci that diverge in allele frequency between extreme temperatures. We've already identified a polymorphic phenotype across this temperature gradient, so it might be interesting to see how well environmentally associated allelic variation lines up with alleles that associate with our phenotype. We can find these phenotypically-associated alleles with a GWA. 


## Fst ##

Since our PCA and structure plots were still a bit mysterious, we also might curious about what Fst among populations shows - for example, when calculating the mean Fst across all pairwise comparisons. 

We are also interested whether this any signatures of selection driven by temperature adaptation, so we might as well do these calculations at a fine-scale - in windows across the genome - to see if we can find any potentially adaptive regions. While there are tons of programs that calclulate Fst (plink, Pixy, SNPrelate), here we're going to use VCFtools.


```bash
mkdir analysis/fst_comparisons

#For vcftools --wier-fst-pop, we need a file containing sample names, for each pop

for pop in `cut -d"." -f2 vcf/Chinook_GWAS_filtered_fixedsamps.fam | uniq`
do
cut -d" " -f1 vcf/Chinook_GWAS_filtered_fixedsamps.fam | grep -w "$pop" > ~/analysis/fst_comparisons/${pop}.samples
done

```

Our input data is all set up for VCFtools but we have to consider how we can get all pairwise comparisons efficiently.
It helps to write out some _psuedocode_.

Comparisons of interest:  
For population 1, 2:10  
For population 2, 3:10  
For population 3, 4:10  
...  
For population 9, 10  

Breaking this down helps outline two major steps 1) that **for each population of interest** (left side), we want to compare **each population that comes after it**. This is a nested loop structure.

```bash
#how can we set our comparison pop to be dependent on the intial focal pop?
#dont forget how to do basic arithmetic in bash!
pop2=`echo $((1+1))`
echo $pop2

#but we want to compare our focal pop not just to the next one, but all pops that follow. seq helps with this.

pop_comp=`seq $((1+1)) 10`
echo $pop_comp

#Now we can put this all together in a nested loop to perform all population comparisons
#testing:
for pop1 in {1..9}
do
        for pop2 in `seq $((pop1+1)) 10`
        do
        echo $pop1 $pop2 
	done
done

```

Now we can use this same logic to run all pairwise population Fst comparisons

```bash
for i in {1..9}
do
	for k in `seq $((i+1)) 10`
	do

	vcftools \
	--gzvcf vcf/Chinook_GWAS_filtered_fixedsamps.vcf.gz \
	--weir-fst-pop ~/analysis/fst_comparisons/p$i.samples \
	--weir-fst-pop ~/analysis/fst_comparisons/p$k.samples \
	--out ~/analysis/fst_comparisons/pop${i}_pop${k}_10kb \
	--fst-window-size 10000 \
	--fst-window-step 10000

	done
done

#list all the log files from the vcftools fst calc runs
ls analysis/fst_comparisons/*.log

#lets check out one
less analysis/fst_comparisons/pop1_pop10_10kb.log
#what we're interested in is the weighted fst estimate
#this accounts for differences in coverage across sites

grep "Weir and Cockerham weighted Fst estimate:" analysis/fst_comparisons/*.log
```

### Code Break Q! ###

Make a results file of the pairwise weighted fst estimates based on the grep command above. Use pipes and basic UNIX commands like _tr_, _cut_,and _sed_ to split the output into a space seperated file with three columns: 1) pop A, 2) pop B, and 3) Fst. Save it as analysis/fst_comparisons/weighted_fst_pairwise.txt

<details>
<summary markdown="span">**Answer**
</summary>
```bash
grep "Weir and Cockerham weighted Fst estimate:" analysis/fst_comparisons/*.log | tr ":" "\t"  | sed 's|analysis/fst_comparisons/||g' | sed 's|_10kb.log||g' | cut -d$'\t' -f1,3 | tr "_" "\t" > analysis/fst_comparisons/weighted_fst_pairwise.txt
```
</details>

If you're not sure how that worked out, break down each part of the pipe, step by step and see how it changes as a result of each new command


## Genome-wide Association ##

A nice thing about plink, the program we used in the previous tutorial, is that alot of programs take the .bed/.fam format for input, including the GWA program GEMMA. We're going to use the same VCF we used to infer patterns of population structure. It's amazing how easy (and fast!!) it is to run a GWA, but we have to be careful about the statistical design and interpretation of these type of analyses.

By default, GEMMA knows to take the 6th column of the plink .fam file as the dependent variable. Take a look at our fam file right now using 'less -S vcf/Chinook_GWAS_filtered_fixedsamps.fam'. Darn, it doesn't have any phenotype info! Let's replace the current column with phenotype info in phenos.txt

```bash
#to do this, lets first make sure we the phenotype information we have is in the same order as our samples
comm --help #make sure you know how to interpret the output!
comm  <( cut -d" " -f1 vcf/Chinook_GWAS_filtered_fixedsamps.fam  ) <( cut -d"," -f-1 /mnt/data/vcf/phenos.txt ) #notice only column 2 (rows in common between datasets) prints
#note the subprocessing

#make modified fam file, by pasting the first 5 columns of our fam with the 2nd column of the phenotype file.
#this will make plink happy as now our 6 column = phenotype information

paste -d " "  <( cut -d" " -f1-5 vcf/Chinook_GWAS_filtered_fixedsamps.fam) <( cut -d"," -f2 /mnt/data/vcf/phenos.txt) > vcf/Chinook_GWAS_filtered_fixedsamps.fammod

```
For operations like above, we're lucky that our information we're interested in merging is in the same order (alphanumeric sample order). Otherwise we'd want to incorporate a couple sort commands in there. This is a super important thing to make sure of before doing any merging

```bash
#plink expects the phenotype to be in the -bfile <prefix>.fam, but right now its in <prefix>.fammod. lets do some quick renaming
mv vcf/Chinook_GWAS_filtered_fixedsamps.fam vcf/Chinook_GWAS_filtered_fixedsamps.famnophenos
mv vcf/Chinook_GWAS_filtered_fixedsamps.fammod vcf/Chinook_GWAS_filtered_fixedsamps.fam

#actually running the GWA is as simple as:
/mnt/software/gemma-0.98.1-linux-static \
	-bfile vcf/Chinook_GWAS_filtered_fixedsamps \
	-lm \
	-o Chinook_GWAS \
	-miss .10 \
	-maf 0.01

#gemme provides options for specifying whether we want to exclude sites based on missing data (-miss) or minor allele frequency (-maf).
```

We talked aobut in lecture how population structure can have effects on GWA that are very important to control for. Lets compare the above GWA to one that controls for population structure via the relatedness matrix in a linear mixed model framework

```r
 #first use the option -gk to generate the relatedness matrix 
/mnt/software/gemma-0.98.1-linux-static \
	-bfile vcf/Chinook_GWAS_filtered_fixedsamps \
	-gk \
	-o Chinook_GWAS_filtered_fixedsamps

#run the GWA controlling for relatedness, with the output matrix (*.cXX.txt)
/mnt/software/gemma-0.98.1-linux-static \
	-bfile vcf/Chinook_GWAS_filtered_fixedsamps \
	-k output/Chinook_GWAS_filtered_fixedsamps.cXX.txt \
	-lmm 4 \
	-o Chinook_GWAS_relatedness
```

GEMMA writes results from these analyses to a folder it makes called output/. Lets rename this folder and move it inside analysis/, and because we are impatient we'll just take a quick look at the the min p-values from these two different GWA runs, before we read it into R.

```bash

mv output/ analysis/gwas/
head analysis/gwas/Chinook_GWAS.assoc.txt #lets look at the p_wald values, column 11
head analysis/gwas/Chinook_GWAS_relatedness.assoc.txt #note that the p_wald column for the linear mixed effect GWA is 13

#a quick way to figure out what SNP is most signficant would just be resorting by those columns, taking just the first row (and dropping all the other columns we don't want to see)

sort -g -k11,11 analysis/gwas/Chinook_GWAS.assoc.txt | head -n 2 | cut -d$'\t' -f11  #g tells sort to interpret scientific notation, -k to sort on column 11
sort -g -k13,13 analysis/gwas/Chinook_GWAS_relatedness.assoc.txt | head -n 2 | cut -d$'\t' -f13

#Note that our minimum p-value is much lower when we account for relatedness

```

Now that we've run some analyses to investigate potential patterns of selection, lets move to R to do some statistics and data visualization.

