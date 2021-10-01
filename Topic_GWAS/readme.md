---
title: "Topics 8: Population genomics - Selection"
permalink: /Topic_GWAS/
topickey: 9
topictitle: "Population genomics - Seletion"
---

## Accompanying material

## Daily Questions:
1. For a site that is invariant in both populations, what is Fst?
2. What is the average Fst between the most distant populations (pop 1 vs pop 10)? What about nearby populations (pop 5 v pop 6?)
3. What effect does uncontrolled population structure often have on GWA qqplot?

Last topic we used bcftools to filter snps, and plink & admixture to investigate population structure. We will plot these results to get a better visual sense of this data, but not quite yet! Before we move to R, we are going to run some analyses to infer potential patterns of selection and adaptation in our dataset. 

Remember, we want to test hypothesis about how our populations of Chinook salmon might have adapted to different temperature environments. We could do this a bunch of different ways, but lets try two. The approach we'll take is to see if we can identify loci that show particularly extreme differentiation between environments of varying temperatures. We'll use the Fst statistic to find these loci that diverge in allele frequency between extreme temperatures. We've already identified a polymorphic phenotype across this temperature gradient, so it might be interesting to see how well environmentally associated allelic variation lines up with alleles that associate with our phenotype. We can find these phenotypically-associated alleles with a GWA. 

Lets start with Fst. To keep things simple for now, we're going to calculate Fst between population 1 and population 10, a comparison that shows the largest difference in temperature. We are going to calculate Fst between these groups using the perl tool vcf2fst.pl. This is a custom script from Greg Owens, since we want to keep the numerator and denominator from the Fst calculations, which is hard to do. (Lots of other programs calculate Fst - can you think of a reason why might it be important to keep track of both the denominator and numerator?)

We need two files, a sample info file and a group file. The sample info file tells the program which population each sample is in and the group file tells the program which populations to compare. We can make them here:

```bash
for i in `cat samplelist.txt`;
        do echo -e "$i\t${i/%????/}";
done > sampleinfo.txt
echo -e "ANN\t1\nARG\t2" > popinfo.txt
```

Next, less run a GWAS. A nice thing about plink is that alot of programs take the .bed/.fam format for input, including the GWA program GEMMA. We're going to use the same VCF we used to infer patterns of population structure. It's amazing how easy it is to run a GWA, but we have to be careful about the statistical design and interpretation of these type of analyses.

By default, gemma knows to take the 6th column of the plink .fam file as the dependent variable. So first, we need to modify this fam file to include or phenotype of interest.

```bash
#modify the fam file, replacing the 6th file with our actual phenotypes
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

GEMMA writes results from these analyses to a folder it makes called 'output'. Because we are impatient we'll just take a quick look at the the min p-values from these two different GWA runs, before we read it into R.

``bash
cd output
head Chinook_GWAS.assoc.txt #lets look at the p_wald values, column 11
head Chinook_GWAS_relatedness.assoc.txt #note that the p_wald column for the linear mixed effect GWA is 13

sort -g -k11,11 Chinook_GWAS.assoc.txt | head -n 2 | cut -d$'\t' -f11  #g tells sort to interpret scientific notation
sort -g -k14,14 Chinook_GWAS_relatedness.assoc.txt | head -n 2 | cut -d$'\t' -f14

#our minimum p-value is much lower when we account for relatedness!
```

We've now run several analyses to investigate population structure and potentially, local adaptation. Lets move to R to do some statistics and data visualization.

