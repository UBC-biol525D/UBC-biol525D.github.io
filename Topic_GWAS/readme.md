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

```bash
#by default, gemma knows to take the 6th column of the plink .fam file as the dependent variable
#so first, we need to modify this fam file to include or phenotype of interest
paste <( cut -d" " -f1-5 vcf/Chinook_GWAS_fiiltered_fixedsamps.fam) phenos.txt

#running the GWA is as simple as:
~/software/gemma-0.98.1-linux-static -bfile vcf/Chinook_GWAS_fiiltered_fixedsamps -lmm 4 -o Chinook_GWAS

#this will take a little bit of time, and write output to /output file of your current working directory
#but remember the effects that uncontrolled population structure can have on GWA type analyses.

#lets compare to the above GWA to one that controls for population structure via the relatedness matrix
~/software/gemma-0.98.1-linux-static -bfile vcf/Chinook_GWAS_fiiltered_fixedsamps -gk -o Chinook_GWAS_fiiltered_fixedsamps #gk is the option for generating the relatedness matrix

#run the GWA controlling for relatedness
~/software/gemma-0.98.1-linux-static -bfile vcf/Chinook_GWAS_fiiltered_fixedsamps -k output/Chinook_GWAS.cXX -lmm 4 -o analysis/Chinook_GWAS_relatedness
```


We've now run several analyses to investigate population structure and potentially, local adaptation. Lets move to R to do some statistics and data visualization.

