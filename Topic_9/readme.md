---
title: "Topics 9: Population genomics - Selection"
permalink: /Topic_GWAS/
topickey: 9
topictitle: "Population genomics - Selection"
---

## Daily Questions:
1. For a site that is invariant in both populations, what is Fst?
2. What is the average Fst between the most distant populations (pop 1 vs pop 10)? What about nearby populations (pop 5 v pop 6?)
3. What effect does uncontrolled population structure often have on GWA qqplot?

Last topic we used bcftools to filter snps, and plink & admixture to investigate population structure. We will plot these results to get a better visual sense of this data, but not quite yet! Before we move to R, we are going to run some analyses to infer potential patterns of selection and adaptation in our dataset. 

Remember, we want to test hypothesis about how our populations of Chinook salmon might have adapted to different temperature environments. We could do this a bunch of different ways, but lets try two. The approach we'll take is to see if we can identify loci that show particularly extreme differentiation between environments of varying temperatures. We'll use the Fst statistic to find these loci that diverge in allele frequency between extreme temperatures. We've already identified a polymorphic phenotype across this temperature gradient, so it might be interesting to see how well environmentally associated allelic variation lines up with alleles that associate with our phenotype. We can find these phenotypically-associated alleles with a GWA. 


## Fst ##

Since our PCA and structure plots were still a bit mysterious, we also might curious about what Fst among populations shows - for example, when calculating the mean Fst across all pairwise comparisons. 

We are also interested whether this any signatures of selection driven by temperature adaptation, so we might as well do these calculations at a fine-scale - in windows across the genome - to see if we can find any potentially adaptive regions. While there are tons of programs that calclulate Fst (plink, Pixy, SNPrelate), here we're going to use VCFtools.


```bash
#For vcftools --wier-fst-pop, we need a file containing sample names, for each pop

for pop in `cut -d"." -f2 vcf/Chinook_GWAS_fiiltered_fixedsamps.fam | uniq`
do
cut -d" " -f1 vcf/Chinook_GWAS_fiiltered_fixedsamps.fam | grep -w "$pop" > ${pop}.samples
done

#Use vcftools to calcluate fst in 10kb windows, across all pop pairs
mkdir analysis/fst_comparisons
#vcftools can't read bgzipped files, but can read gzipped
gunzip vcf/Chinook_GWAS.filtered.fixedsamps.vcf.gz vcf/Chinook_GWAS.filtered.fixedsamps.vcf
gzip vcf/Chinook_GWAS.filtered.fixedsamps.vcf

for i in {1..9}
do
	for k in `seq $((i+1)) 10` # $((EXPR)) is bash for arithmetic expression!
	do
	vcftools --gzvcf vcf/Chinook_GWAS.filtered.fixedsamps.vcf.gz --weir-fst-pop p$i.samples --weir-fst-pop p$k.samples --out analysis/fst_comparisons/pop${i}_pop${k}_10kb --fst-window-size 10000 --fst-window-step 10000
	done
done

ls analysis/fst_comparisons/*.log
grep "Weir and Cockerham weighted Fst estimate:" analysis/fst_comparisons/*.log
```

### Code Break Q! ###

Make a results file of the pairwise weighted fst estimates based on the grep command above. Use pipes and basic UNIX commands like _tr_, _cut_,and _sed_ to split the output into a space seperated file with three columns: 1) pop A, 2) pop B, and 3) Fst. Save it as analysis/fst_comparisons/weighted_fst_pairwise.txt

SOLUTION (don't click me unless your stuck!):
```bash
grep "Weir and Cockerham weighted Fst estimate:" analysis/fst_comparisons/*.log | tr ":" "\t"  | sed 's|analysis/fst_comparisons/||g' | sed 's|_10kb.log||g' | cut -d$'\t' -f1,3 | tr "_" "\t" > analysis/fst_comparisons/weighted_fst_pairwise.txt
```
{: .spoiler}



## Genome-wide Association ##

A nice thing about plink, the program we used in the previous tutorial, is that alot of programs take the .bed/.fam format for input, including the GWA program GEMMA. We're going to use the same VCF we used to infer patterns of population structure. It's amazing how easy it is to run a GWA, but we have to be careful about the statistical design and interpretation of these type of analyses.

By default, GEMMA knows to take the 6th column of the plink .fam file as the dependent variable. So first, we need to modify this fam file to include or phenotype of interest.

```bash
#modify the fam file, replacing the 6th column with our actual phenotypes
#first, check whether they are in the same order using comm (how do we interpret the output?)
comm  <( cut -d" " -f1 vcf/Chinook_GWAS_fiiltered_fixedsamps.fam | sort ) <( cut -d"," -f-1 phenos.txt | sort)

#make modified fam file
paste -d " "  <( cut -d" " -f1-5 vcf/Chinook_GWAS_fiiltered_fixedsamps.fam) <( cut -d"," -f2 phenos.txt) > vcf/Chinook_GWAS_fiiltered_fixedsamps.fammod

#plink expects the phenotype to be in the -bfile <prefix>.fam, but right now its in <prefix>.fammod. lets do some quick renaming
mv vcf/Chinook_GWAS_fiiltered_fixedsamps.fam vcf/Chinook_GWAS_fiiltered_fixedsamps.famnophenos
mv vcf/Chinook_GWAS_fiiltered_fixedsamps.fammod vcf/Chinook_GWAS_fiiltered_fixedsamps.fam

#actually running the GWA is as simple as:
~/software/gemma-0.98.1-linux-static -bfile vcf/Chinook_GWAS_fiiltered_fixedsamps -lm -o Chinook_GWAS -miss .10 -maf 0.01

#lets compare the above GWA to one that controls for population structure via the relatedness matrix in a linear mixed model framework
~/software/gemma-0.98.1-linux-static -bfile vcf/Chinook_GWAS_fiiltered_fixedsamps -gk -o Chinook_GWAS_fiiltered_fixedsamps #gk is the option for generating the relatedness matrix
#run the GWA controlling for relatedness
~/software/gemma-0.98.1-linux-static -bfile vcf/Chinook_GWAS_fiiltered_fixedsamps -k output/Chinook_GWAS_fiiltered_fixedsamps.cXX.txt -lmm 4 -o Chinook_GWAS_relatedness
```

GEMMA writes results from these analyses to a folder it makes called output/. Lets rename this folder and move it inside analysis/, and because we are impatient we'll just take a quick look at the the min p-values from these two different GWA runs, before we read it into R.

```bash

mv output/ analysis/gwas/
head analysis/gwas/Chinook_GWAS.assoc.txt #lets look at the p_wald values, column 11
head analysis/gwas/Chinook_GWAS_relatedness.assoc.txt #note that the p_wald column for the linear mixed effect GWA is 13

sort -g -k11,11 analysis/gwas/Chinook_GWAS.assoc.txt | head -n 2 | cut -d$'\t' -f11  #g tells sort to interpret scientific notation
sort -g -k14,14 analysis/gwas/Chinook_GWAS_relatedness.assoc.txt | head -n 2 | cut -d$'\t' -f14

#our minimum p-value is much lower when we account for relatedness..

```

Now that we've run some analyses to investigate potential patterns of selection, lets move to R to do some statistics and data visualization.

