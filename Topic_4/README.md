---
title: "Topic 3n4: Genome Assembly"
permalink: /Topic_4/
topickey: 4
topictitle: "Genome Assembly"
---

## Recorded lecture

<iframe src="https://monash.au.panopto.com/Panopto/Pages/Embed.aspx?id=835dc4b6-e73c-4892-88cc-ac820189a883&autoplay=false&offerviewer=true&showtitle=true&showbrand=false&start=0&interactivity=all" height="405" width="720" style="border: 1px solid #464646;" allowfullscreen allow="autoplay"></iframe>

## Accompanying material
* Slides 2020: [UBC - De novo Assembly 2020](./Assembly2020.pdf)
* Background reading: The present and future of de novo whole-genome assembly [Paper](https://academic.oup.com/bib/article/19/1/23/2339783?login=true#119542667). 

Programming Resources
* [Flye Manual](https://github.com/fenderglass/Flye/blob/flye/docs/USAGE.md)
* Flye Paper. Assembly of Long Error-Prone Reads Using Repeat Graphs [Paper](https://doi.org/10.1038/s41587-019-0072-8).
* [Hslr Manual](https://github.com/vpc-ccg/haslr)
* Hslr Paper. HASLR: Fast Hybrid Assembly of Long Reads [Paper](https://www.cell.com/iscience/fulltext/S2589-0042(20)30577-0?_returnURL=https%3A%2F%2Flinkinghub.elsevier.com%2Fretrieve%2Fpii%2FS2589004220305770%3Fshowall%3Dtrue).
* [Spades Manual](https://github.com/fenderglass/Flye/blob/flye/docs/USAGE.md)
* Spades Paper. A New Genome Assembly Algorithm and Its Applications to Single-Cell Sequencing [Paper](https://cab.spbu.ru/files/release3.15.3/manual.html).


# Topic 3 Genome Assembly

# Code break questions


1) How long are the sequences in /home/biol525d/data/shortreads/SalmonSim.Stabilising.p1.1.6400000_R1.fastq.g? what about in /home/biol525d/data/longreads/SalmonSim.Stabilising.p1.3.30k.PacBio.fastq.gz?

> Hint: One approach could be to subset the fastq file to just retain lines with the actually sequence info. Play around with these commands.
```bash
cat SalmonSim.Stabilising.p1.1.6400000_R1.fastq.gz | head #why doesn't this work?
zcat SalmonSim.Stabilising.p1.1.6400000_R1.fastq.gz | grep "chr" -A1 #what does the A option do?
#what could you add to this pipe to only keep the sequence length? hint: egrep -v "match1|match2"
#pipe the final output that is just the read sequence lines to the following awk command. 
awk -F "" '{print NR "\t" NF}'  #what is NR and what is NF? what does the -F "" option do?
```
2) Given the genome size of 10Mb that we will be working with, what is the expected estimate of coverage/bp in just the forward reads?
> Hint: to get mean of column two you can try:  `awk '{ total += $2 } END { print total/NR }' ` \
> Hint:  `wc -l` gives a count of the number of lines \
> Note: math is kind of annoying in command lind. expr works pretty well though for basic purposes (make sure to leave spaces around operators).
> you can also do cool things like save the output of commands to a variable using "` `"
```bash
readlength=`zcat SalmonSim.Stabilising.p1.1.6400000_R1.fastq.gz | grep "chr" -A1 |  egrep -v "chr|--" |\
 awk -F "" '{print NR "\t" NF}' | awk '{ total += $2 } END { print total/NR }' `
numreads=`zcat SalmonSim.Stabilising.p1.1.6400000_R1.fastq.gz | grep "chr" -A1 |  egrep -v "chr|--" | wc -l`
echo "mean coverage =" "$(( ($readlength * $numreads) / 10000000 ))"
```

## Tutorial 

Your goal for today is to see how well we can assemble a portion of the Salmon reference genome with short read, short+long read (hybrid), and long read only approaches. New software with improved algorithms and sequence technology-specific-approaches are constantly being developed so when you start your own assembly project, make sure to look around for the newest and most promising methods. It will also be highly dependent on the architecture of the genome you are trying to assemble (i.e. depending on repeat content, ploidy, etc.). Beyond just running the commands to run an assembler, we're going to focus on the things you should be considering leading up to assembling your genome (preproccessing/data quality assesment) and after (how good is my assembly?!).

You've already spent some time getting familiar with the data we're working with - Our overall aim is to understand the evolutionary processes that allow salmon populations to persist across a key temperature gradient, and having a high quality reference genome will be indispensible for looking at not just the the frequency of variant nucleotides but their position relative to one another and other features of the genome. 


We already used some simple command line tools to look at read length, number of reads, and expected coverage above. However, like most bioinformatics, there are a bunch of tools that have been developled to look at various statistics related to your raw reads and we don't need to reinvent the wheel. We are going to use *fastqc* to get a better sense of the quality of our read data.

```bash
mkdir fastqc && cd fastqc
```

For our assembly, since we only sequence one individual, we have two files for our short reads (one for forward reads and one for reverse) and one for long reads. We're going to be referring to these files alot so lets first make a short cut path to these files and then run fastqc on each.

## Check the quality of your reads

```bash
shortreads="/home/biol525d/data/shortreads/"
longreads="/home/biol525d/data/longreads/"

fastqc ${shortreads}/*fastq.gz ${longreads}/*fastq.gz -o ./ 
#some programs allow you to input multiple files at once using wild cards. this is convienient. 
#otherwise we would have had to specify each file individually. heres anotherway we could use wildcards here: 
#fastqc ${shortreads}/*R?.fastq.gz ${longreads}/*fastq.gz -o ./ )

#summarise into a single report
multiqc ./
```

Download the output of multiqc run as follows:
```bash
scp <username@ip.address>:~/Topic_3n4/fastqc/multiqc_report.html <path on your computer where you want the file>
```

Since we are using simulated data, our data does not have adapters. But taking a look at the signatures of our highly quality reads can help you get a sense of what very high quality data might look like, and additionally,by comparing to your real data, how sequencing technology influences data quality.  

Discussion Quesiton: our data effectively has had adapters removed from reads and already been trimmed for low quality BPs near the end of reads. what might be the cons of read trimming? 

## Let's go ahead with genome assembly

Genome assembly can take a long time. Because our course is short we're only going to actually generate an assembly from the program that runs fastest, but we're going to provide the commands for several programs that we've already run, and together we will compare their outputs. 


#### Short read assembly: SPADES - *don't run*
#using a single core this takes ~14 minutes
```bash
#mkdir spades
#/ohta1/julia.kreiner/software/SPAdes-3.
#15.3-Linux/bin/spades.py \
	--pe1-1 ${shortreads}/SalmonSim.Stabilising.p1.1.6400000_R1.fastq.gz \
	--pe1-2 ${shortreads}/SalmonSim.Stabilising.p1.1.6400000_R2.fastq.gz \
	-o ./spades/
````


#### Hybrid (short + long) assembly - SPADES & HASLR - *don't run*
#haslr takes ~15 minutes
```bash
#make new dir
mkdir hybridspades
#run
#/ohta1/julia.kreiner/software/SPAdes-3.15.3-Linux/bin/spades.py \
	--pe1-1 ${shortreads}/SalmonSim.Stabilising.p1.1.6400000_R1.fastq.gz \
	--pe1-2 ${shortreads}/SalmonSim.Stabilising.p1.1.6400000_R2.fastq.gz \
	--pacbio ${longreads}/SalmonSim.Stabilising.p1.3.30k.PacBio.fastq.gz \
	-t 20 \
	-o /hybridspades/ 

#mkdir hybridhaslr
#install
#conda install -c bioconda haslr
#run
#haslr.py \
	-t 20 \
	-s ${shortreads}/SalmonSim.Stabilising.p1.1.6400000_R1.fastq.gz ${shortreads}/SalmonSim.Stabilising.p1.1.6400000_R2.fastq.gz \
	-l ${longreads}/SalmonSim.Stabilising.p1.3.30k.PacBio.fastq.gz \
	-x pacbio \
	-g 10m \
	-o hybridhaslr \
	--cov-lr 40 #takes 15 minutes
#
```

#### Long read assembly - FLYE - *don't run*
#this took 8 minutes with 20 threads

```bash
#install
#conda install flye

#run flye assuming lower quality pacbio reads
flye \
	--pacbio-raw ${longreads}/SalmonSim.Stabilising.p1.3.30k.PacBio.fastq.gz \
	--threads 20 \
	-o flye/ 
	--genome-size 10m
```

Now we have four assemblies, 1 short read (spades), 2 hybrid (spades, haslr), and one long-read (flye).

## Assess quality of assemblies

Lets compare assemblies.

```bash
mv spades/scaffolds.fasta spades.fasta
mv hybridspades/scaffolds.fasta spades_hybrid.fasta
mv hybridhaslr/asm_contigs_k49_a3_lr40x_b500_s3_sim0.85/asm.final.fa haslr_hybrid.fasta
mv flye/assembly.fasta flye_longread.fasta
```

bbmap is command line alignment program that has a collection of nice scripts for library and sequence quality control. We're going to use its stats.sh script to get at some basic stats related to the number and length of sequences in our assembly.

An important stat is N50: the number of contigs making up 50% of your assembly. 
Relatedly, L50 describes for those largest sequences that make up 50% of your assembly, what the minimum sequence length is.
This becomes more intutitive when you look at the lower table outputted by bbmap (unfortunately there is some inconsitency between these definitions in that L/N usage is often swapped, as you'll see later)

We have several assemblies we want to generate stats for, so lets automate this using a for loop:
```bash
for i in haslr flye haslr hybridspades spades 
do

~/software/bbmap/stats.sh in=${i}.fa > ${i}.stats

done
```

Question: which genome has the best *contiguity* via N50/L50 metrics? What might be an issue with relying on only contiguity? 

Having a reference genome of a closely related species can really help asses how well we've done with assembly, outside of just contiguity. In our case we have one better - the high quality reference genome from which our reads were simulated. Lets test out quast - a program for getting detail stats when comparing to a closely related high quality reference.

Importantly, this allows us to get not just the number and length of our contigs/scaffolds, but also _completeness_ (possible when you know the target genome size) and _correctness_.

```bash
#it would be nice to know how well our assemblies have captured known genes. Lets extract this information from the true reference genome's gff file and format it as quast expects 
#quast expects the column order as follows: chromosome, gene id, end, start
awk 'split($9,a,";") {print $1 "\t" a[1] "\t" $5 "\t" $4}' ../data/SalmonAnnotations_forIGV.gff | sed 's/ID=//g' > SalmonReference.genes

#run quast on all 4 assemblies at once
mkdir quast
/ohta1/julia.kreiner/software/quast/quast.py flye.fa haslr.fa hybridspades.fa spades.fa \
	-r SalmonReference.fasta \
	-g SalmonReference.genes \
	-o quast

scp -r <username@ip.address>:~/Topic_3n4/quast/ ./
```

Open the report.html results file in your browser and explore the outcome of our assembly efforts. Make sure to eventually click the "View in icarus contig browser".

Question: what new information does this report give us? how does it inform on completeness and correctness?



## Review some important command line operations we used today.

For loops - this can be a list of items seperated by spaces or you can specify of range numeric values to use a an iterator
```bash
for i in {1..100}
do
echo $i
done
```

Making new folders `mkdir` \
Renaming files `mv` and moving them  `mv file path/new_file_name` \
Counting `wc -l` \
Find and replace `sed 's/find/replace/g'` \
Column means with `awk '{ total += $2 } END { print total/NR }'` \
Assigning variables `shortreads="/home/biol525d/data/shortreads/"` and calling them `echo ${shortreads}` \

Printing and splitting certain columns (specified with $) with awk 
```bash
awk 'split($9,a,";") {print $1 "\t" a[1] "\t" $5 "\t" $4}'
#split column 9 by the delimeter ";" and save every split value into the array a
#print column 1, the first value of array a, column 5, and column 4
```

Cut also allows your to print certain columns but the order will be the same as provided in the input. For e.g., compare the following outputs.
```bash
cut -f2,1 SalmonReference.genes 
cut -f1,2 SalmonReference.genes
```
