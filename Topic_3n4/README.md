---
title: "Topic 3n4: Genome Assembly"
permalink: /Topic_3n4/
topickey: 5
topictitle: "Genome Assembly"
---

## Recorded lecture

<iframe src="https://monash.au.panopto.com/Panopto/Pages/Embed.aspx?id=835dc4b6-e73c-4892-88cc-ac820189a883&autoplay=false&offerviewer=true&showtitle=true&showbrand=false&start=0&interactivity=all" height="405" width="720" style="border: 1px solid #464646;" allowfullscreen allow="autoplay"></iframe>

## Accompanying material

* Slides 2020: [UBC - De novo Assembly 2020](./Assembly2020.pdf)
* Background reading: Comparison of the two major classes of assembly algorithms: overlap-layout-consensus and de-bruijn-graph [Paper](./background_reading/Briefings in Functional Genomics-2011-Li-bfgp-elr035.pdf). Briefings in Functional Genetics 2011.
* Sense from sequence reads: methods for alignment and assembly [Paper](./background_reading/Flicek&Birney2009.pdf). Flicek and Birney. Nature methods supplement 2009.
* [Flye Manual](https://github.com/fenderglass/Flye/blob/flye/docs/USAGE.md)
* Flye Paper. Assembly of Long Error-Prone Reads Using Repeat Graphs [Paper](https://doi.org/10.1038/s41587-019-0072-8).


# Topic 5 Genome Assembly

# Code break questions


1) how long are the sequences in /home/biol525d/data/shortreads/SalmonSim.Stabilising.p1.1.6400000_R1.fastq.g? what about in /home/biol525d/data/longreads/SalmonSim.Stabilising.p1.3.30k.PacBio.fastq.gz?

Hint: One approach could be to subset the fastq file to just retain lines with the actually sequence info. Test out these commands.
```bash

cat SalmonSim.Stabilising.p1.1.6400000_R1.fastq.gz | head #why doesn't this work?
zcat SalmonSim.Stabilising.p1.1.6400000_R1.fastq.gz | grep "chr" -A1 #what does the A option do?
#what could you add to this pipe to only keep the sequence length? hint: egrep -v "match1|match2"
#pipe the final output that is just the read sequence lines to the following awk command. 
awk '{print NR "\t" NF}'  #what is NR and what is NF?

```
2) Given a genome size of 10Mb (10,000,000 bps), what do is your estimate of coverage/bp in just the forward reads?

Hint: to get mean of column two you can try:  `awk '{ total += $2 } END { print total/NR }' `
Hint:  `wc -l` gives a count of the number of lines

Note: math is kind of annoying in command lind. expr works pretty well though for basic purposes (make sure to leave spaces around operators).

## Tutorial 

Your goal for today is to see how well we can assemble a portion of the Salmon reference genome with short read, short+long read, and long read only approaches. New software with improved algorithms and sequence technology specific approaches so we might be careful about generalizing "the best assembly approach" but it is worthwile getting familiar with some of these.

You've already spent some time getting familiar with the data we're working with - Salmon populations that have been evolving to a temperature gradient - however if we want to learn about the evolutionary processes going on in these pops its going to be really useful to have a genome to align it to, to get information about not just variant nucleotides but their relative positioning. 


Details for this dataset are as follows: 

Species:

Actual genome size: 10Mb

Type: Illumina Paired End, Pacbio

Read number (total - including both reads per pair):


Read size (each read):

Insert length (sd): 

We used some simple command line tools to look at read length above. However, like most bioinformatics, there are a bunch of tools that have been developled to look at various statistics related to your raw reads. We are going to use fastqc to get a better sense of our data.

```bash
mkdir fastqc && cd fastqc
```

For our assembly, since we only sequence one individual, we have two files for our short reads (one for forward reads and one for reverse) and one for long reads. We're going to be referring to these files alot so lets first make a short cut path to these files and then run fastqc on each.

## Check the quality of your reads

```bash
shortreads="/home/biol525d/data/shortreads/"
longreads="/home/biol525d/data/longreads/"

fastqc ${shortreads}/*fastq.gz ${longreads}/*fastq.gz -o ./ #some programs allow you to input multiple files at once using wild cards. this is convienient. otherwise we would have had to specify each file individually. heres anotherway we could have specified. #Q: what does the option  "-o ./" do?

fastqc ${shortreads}/*R?.fastq.gz ${longreads}/*fastq.gz -o ./ #Q: what does the "?" do? what about the "*"?

#summarise into a single report
multiqc ./
```

Download these the output of these files as follows:
```bash
scp <username@ip.address>:~/Topic_3/fastqc/multiqc_report.html <path on your computer where you want the file>
```

*** I'm not sure how useful this portion will really be - might be worth just discussing that the adapters have already been removed and we have pretty high quality reads, even at the ends ****

Discussion Question: what are the main benefits and weeknesses of the different data types (i.e. pacbio versus illumina short read)
Discussion Quesiton: there are tradeoffs of trimming short read data, to remove low quality bases near the end of the read. 

## Go ahead with genome assembly

Genome assembly can take a long time. Because our course is short we're only going to actually generate an assembly from a program that runs quickly, but we're going to provide the commands for several programs and compare their outputs, which we have provided. 


#short read assembly: SPADES
#using a single core this takes ~14 minutes
```bash
mkdir spades
/ohta1/julia.kreiner/software/SPAdes-3.
15.3-Linux/bin/spades.py --pe1-1 ${shortreads}/SalmonSim.Stabilising.p1.1.6400000_R1.fastq.gz --pe1-2 ${shortreads}/SalmonSim.Stabilising.p1.1.6400000_R2.fastq.gz -o ./spades/
````


#two types of hybrid (short + long) assembly: SPADES & HASLR
```bash
#make new dir
mkdir hybridspades
#run
/ohta1/julia.kreiner/software/SPAdes-3.15.3-Linux/bin/spades.py --pe1-1 ${shortreads}/SalmonSim.Stabilising.p1.1.6400000_R1.fastq.gz --pe1-2 ${shortreads}/SalmonSim.Stabilising.p1.1.6400000_R2.fastq.gz --pacbio ${longreads}/SalmonSim.Stabilising.p1.3.30k.PacBio.fastq.gz -t 20 -o /hybridspades/ 

mkdir hybridhaslr
#install
conda install -c bioconda haslr
#run
haslr.py -t 20 -s ${shortreads}/SalmonSim.Stabilising.p1.1.6400000_R1.fastq.gz ${shortreads}/SalmonSim.Stabilising.p1.1.6400000_R2.fastq.gz -l ${longreads}/SalmonSim.Stabilising.p1.3.30k.PacBio.fastq.gz -x pacbio -g 10m -o hybridhaslr --cov-lr 40 #takes 15 minutes
#
```

#long read assembly - using FLYE - *dont run*

```bash
#install
conda install flye

#run flye assuming lower quality pacbio reads
flye --pacbio-raw  ${longreads}/SalmonSim.Stabilising.p1.3.30k.PacBio.fastq.gz --threads 20 -o flye/ --genome-size 10m
#this took 8 minutes with 20 threads
```

Now we have four assembles, 1 short read (spades), 2 hybrid (spades, haslr), and one long read only (flye). Lets compare assemblies.

```bash
mv spades/scaffolds.fasta spades.fasta
mv hybridspades/scaffolds.fasta spades_hybrid.fasta
mv hybridhaslr/asm_contigs_k49_a3_lr40x_b500_s3_sim0.85/asm.final.fa haslr_hybrid.fasta
mv flye/assembly.fasta flye_longread.fasta
```

bbmap is command line alignment program that has a collection of nice scripts for library and sequence quality control. We're going to use its stats.sh script to get at some basic stats around the number and length of sequences in our assembly.

An important stat is N50: the number of contigs making up 50% of your assembly. 
Relatedly, L50 describes for those largest sequences that make up 50% of your assembly, what the minimum sequence length is.
This becomes more intutitive when you look at the lower table outputted by bbmap.

We have several assemblies we want to generate stats for, so lets automate this using a for loop
```bash
for i in haslr flye haslr hybridspades spades 
do

~/software/bbmap/stats.sh in=${i}.fa > ${i}.stats

done
```

Question: which genome has the best *contiguity* via N50/L50 metrics? What might be an issue with relying on only contiguity? 

Having a reference genome of a closely related species can really help asses how well we've done, outside of just contiguity. In our case we have one better - the high quality reference genome from which our reads were simulated. Lets test out quast - a program for getting detail stats when comparing to a closely related high quality reference.

Importantly, this allows us to get at both completness (possible when you know the target genome size) and correctness.


```bash
#for our purposes, it will be useful to know the location of known genes. Lets format this as quast expects 
#column order as follows: chromosome, gene id, end, start
awk 'split($9,a,";") {print $1 "\t" a[1] "\t" $5 "\t" $4}' ../data/SalmonAnnotations_forIGV.gff | sed 's/ID=//g' > SalmonReference.genes

#run quast on all 4 assemblies at once
mkdir quast
/ohta1/julia.kreiner/software/quast/quast.py flye.fa haslr.fa hybridspades.fa spades.fa -r SalmonReference.fasta -g SalmonReference.genes -o quast

scp -r <username@ip.address>:~/Topic_3n4/quast/ ./
```

Open the report.html results file in your browser and explore the outcome of our assembly efforts. Make sure to eventually click the "View in icarus contig browser".



