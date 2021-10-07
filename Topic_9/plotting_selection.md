---
title: "Topic 9 (continued: Selection in R)"
topickey: 9.2
topictitle: "Plotting GWAS/Fst in R"
---
  
Just like before, we're going to need download the output files of the GWA and FST analyses from the server to our local computers.

Use scp again to copy over the updated analysis folder that contains the fst_comparisons/ and the gwas/ folders.


Lets start with *Fst*. 

We have two different outputs that we want to visualize: 
  1) fine scale fst in 10kb windows across the genome
  2) mean fst between all pairwise population comparisons (just a bit more insight into pop structure!)
  
  
Lets plot the latter first since we're still curious about whats going on after the structure and PCA plots.

```r
fst_pairwise<- read.table("analysis/fst_comparisons/weighted_fst_pairwise.txt")
names(fst_pairwise)<-c("x","y","Fst")
head(fst_pairwise)

#it would be nice to visualized the values in a matrix, but we'll have to fool around with the format of the data first
fst_pairwise$x<-as.numeric(gsub("pop","",fst_pairwise$x))
fst_pairwise$y<-as.numeric(gsub("pop","",fst_pairwise$y))

#resort on both columns
fst_pairwise_ordered<-fst_pairwise[with(fst_pairwise, order(x, y)), ]
#if you're not sure why the "with" fxn was important, check out what it does below
with(fst_pairwise, order(x, y))

head(fst_pairwise_ordered)
tail(fst_pairwise_ordered)
#for the matrix, it would be nice if there was a 1,1 and a 10,10 comparison, to square things off. lets add them.
fst_pairwise_ordered <- rbind(c(1,1,NA), fst_pairwise_ordered)
fst_pairwise_ordered <- rbind(fst_pairwise_ordered, c(10,10,NA) )


#OK, back to tidyverse to convert to a matrix
?tidyr:::pivot_wider
fst_wide<-pivot_wider(fst_pairwise_ordered, names_from=y, values_from=Fst)
fst_wide

#and finally, remove the row and column that isn't Fst values, and change class = matrix
fst_wide<-fst_wide[,-1]

```

That took a good amount of work to get the pairwise Fst comparisons into matrix format. Now plotting is pretty simple...

```r
install.packages("pheatmap")
library(pheatmap)
pheatmap(fst_wide,cluster_rows=F, cluster_cols=F, na_col="white",main = "Pairwise Fst")
```

This plot shows really clear isolation by distance, with populations sampled far apart (e.g. 1 and 10) showing the most differentiation and populations sampled closer together being more similar. Super cool.

Now we can dig into the finescale patterns of Fst across the genome. Since Pop 1 and 10 are found at the two extremes of the temperature gradient, lets focus on theme(
  

```r
fst_10kb<-read.table("analysis/fst_comparisons/pop1_pop10_10kb.windowed.weir.fst",header=T)
head(fst_10kb)

#CHROM BIN_START BIN_END N_VARIANTS WEIGHTED_FST    MEAN_FST
#1 chr_1         1   10000         35    0.0994970  0.06130720
#2 chr_1     10001   20000         29    0.0693105  0.04144450
#3 chr_1     20001   30000         39    0.1192490  0.04892740
#4 chr_1     30001   40000         16    0.0556863  0.02148930
#5 chr_1     40001   50000         36    0.0486619  0.02972630
#6 chr_1     50001   60000         26   -0.0136650 -0.00660728
```

We can see that VCFtools kept track of how many SNPs were used to calculate Fst in each window. Lets check this distribution out.

```r
hist(fst_10kb$N_VARIANTS)
summary(fst_10kb$N_VARIANTS) #looks good. nothing below 15 which isn't unreasonable diveristy (15/10000 = 0.0015)

```

We won't do any window filtering here, but variable data quality across windows is something to watch out for. Let's go ahead and plot.default(

```r
ggplot(data=fst_10kb, aes(BIN_START, MEAN_FST)) +
  geom_point() +
  geom_path() +
  theme_bw() +
  facet_wrap(~CHROM)

#this looks a bit too spikey, lets try fitting a non-paramatric loess regression line
ggplot(data=fst_10kb, aes(BIN_START, MEAN_FST)) +
  geom_point() +
  geom_smooth(method="loess",span=.05) +
  theme_bw() +
  facet_wrap(~CHROM)

```

While there are some windows with big Fst values (e.g. near 3e06), we're not seeing any super obvious regions popping out from a basic fst scan. Another way we will look at this is through Genome-wide association.


Rather than looking for allele-frequency differences between populations of interest, GWAS allows us to test whether phenotypes across individuals are well explained by genotypes across individuals, at each locus.

*GWAS*

```r
gwas<-read.table("analysis/gwas/Chinook_GWAS.assoc.txt",header=T)
head(gwas)

```

Note the p_wald column. These are our p-values. GWA has big issues with multiple testing, and also confounding with population structure.
A qq-plot (quantile-quantile plot) allows us to investigate our p-value distribution, which is influenced by both of these things potentially inflating false positives.

```r
install.packages("qqman")
library(qqman)

qq(gwas$p_wald)
```

The observed p-values start deviating from expected values, even at smaller -log10 pvalues. This is the issue we were worried about - an overrepresentation of false positives. 
The one nice thing about our qqplot is that there is an uptick in this difference at our highest pvalues, which might indicate we have a few true positives.
Remember, we also ran a GWAS analysis where we explitly controlled for population structure by including the relatedness matrix as a mixed effect variable. We'll read this in now and compare our two qqplots.

```r
gwas_rltdns<-read.table("analysis/gwas/Chinook_GWAS_relatedness.assoc.txt",header=T)
qq(gwas_rltdns$p_wald)

```

Notice how much how much closer our observed values are to our expected values, especially on the left side of the plot. Unfortunatly for us, though, we don't see any sign of potential true positives. This might be because their is no adaptation, because population structure is confonded with our selective signal, or because we do not have power given our sample size.
It is still informative to see where strongly associated loci lie in the genome. Lets do this for our uncorrected gwas where it looks like we might have a couple hits. 

We will control for multiple testing using a false-discovery rate correction that uses a similar logic as the qqplot.

```r
gwas$p_fdr<-p.adjust(gwas$p_wald,method = "fdr")
fdr05_uncor<-max(gwas$p_wald[gwas$p_fdr < 0.05])

p_uncor <- ggplot(data=gwas, aes(ps,-log10(p_wald))) +
  geom_point() +
  geom_hline(yintercept=-log10(fdr05_uncor), lty="dashed") +
  facet_wrap(~chr) +
  theme_bw()

```

We have two hits, one around 3.4 Mb, and another around 4.5 Mb (that has a few neighbouring friends). 
These SNPs fall within gene X and gene Y, that function in ... It seems convincing that these genes may be involved in thermotolerance, but remember that it is largely counfounded by population structure (our qqplot shows no excess of observed p-values).

It might be of interest to compare the Fst plot, and the GWAS plots before and after correcting for relatedness.

Plotting Challenge
-------------------
Plot the Fst and GWA plots together, vertically aligned using ggarrange.

