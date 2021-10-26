---
title: "Topic 4: Genome Assembly"
permalink: /Topic_4/
topickey: 4
topictitle: "Genome Assembly"
---

## Recorded lecture

<iframe src="https://monash.au.panopto.com/Panopto/Pages/Embed.aspx?id=835dc4b6-e73c-4892-88cc-ac820189a883&autoplay=false&offerviewer=true&showtitle=true&showbrand=false&start=0&interactivity=all" height="405" width="720" style="border: 1px solid #464646;" allowfullscreen allow="autoplay"></iframe>

## Accompanying material
* Slides 2020: [UBC - De novo Assembly 2020](./Assembly2021Julia - JK.pdf)
* Background reading: The present and future of de novo whole-genome assembly [Paper](https://academic.oup.com/bib/article/19/1/23/2339783?login=true#119542667). 

Programming Resources
* [Flye Manual](https://github.com/fenderglass/Flye/blob/flye/docs/USAGE.md)
* Flye Paper. Assembly of Long Error-Prone Reads Using Repeat Graphs [Paper](https://doi.org/10.1038/s41587-019-0072-8).
* [Hslr Manual](https://github.com/vpc-ccg/haslr)
* Hslr Paper. HASLR: Fast Hybrid Assembly of Long Reads [Paper](https://www.cell.com/iscience/fulltext/S2589-0042(20)30577-0?_returnURL=https%3A%2F%2Flinkinghub.elsevier.com%2Fretrieve%2Fpii%2FS2589004220305770%3Fshowall%3Dtrue).
* [Spades Manual](https://github.com/fenderglass/Flye/blob/flye/docs/USAGE.md)
* Spades Paper. A New Genome Assembly Algorithm and Its Applications to Single-Cell Sequencing [Paper](https://cab.spbu.ru/files/release3.15.3/manual.html).


## Code break questions

1. Write a one liner to find all the overlaps (i.e. the beginning or the end) exactly 4bp in length between CTCTAGGCC and a list of other sequences in the file /mnt/data/codebreaks/overlaps.fa

2. Find all the unique 9mers in a fasta sequence /mnt/data/codebreaks/kmer.fa

This might be a tricky one and there are likely many ways to do this. First try out these commands that might help you complete this task. It might also be helpful to determine how many characters are in the sequence (wc -c).

Hints: test out the following commands:

```bash
cut -c1- /mnt/data/codebreaks/kmer.fa
cut -c1-4 /mnt/data/codebreaks/kmer.fa
```

```bash
wc -c /mnt/data/codebreaks/kmer.fa
```

```bash
for num in {1..10}
    do
        echo $num >> file.txt
    done
```
What do these commands do? Can you use commands like this to find all the kmers in the sequence? 


One approach we could use involves variable assignment. Variable assignment is a super powerful part of bash coding. A variable can be used as a shortcut for long path, or it can even contain the output of a command. The notation is as follows:

```bash
shortpath="/mnt/data/shortreads/"
ls $shortpath
cmdout=`echo "test"`
echo $cmdout
```

Think about how you could incorporate some basic algebra and variable assignment to solve this problem.

```bash
#one way to do math in bash is by piping to bc (basic calculator).
echo "1+1" | bc
#another way to do arithmitic in bash is through the $(( )) syntax which tells shell to evaluate its contents
echo $((1+1))
```

<details>
<summary markdown="span">**Answer**
</summary>
```bash
for i in {1..52} 
	do k=$(($i+8))
	cut -c $i-$k /mnt/data/codebreaks/kmer.fa
done
```
</details>


## Tutorial 

Your goal for today is to see how well we can assemble a portion of the Salmon reference genome with short read, short+long read (hybrid), and long read only approaches. New software with improved algorithms and sequence technology-specific-approaches are constantly being developed so when you start your own assembly project, make sure to look around for the newest and most promising methods. It will also be highly dependent on the architecture of the genome you are trying to assemble (i.e. depending on repeat content, ploidy, etc.). Beyond just running the commands to run an assembler, we're going to focus on the things you should be considering leading up to assembling your genome (preproccessing/data quality assesment) and after (how good is my assembly?!).

You've already spent some time getting familiar with the data we're working with - Our overall aim is to understand the evolutionary processes that allow salmon populations to persist across a key temperature gradient, and having a high quality reference genome will be indispensible for looking at not just the the frequency of variant nucleotides but their position relative to one another and other features of the genome. 

For our assembly, since we only sequence one individual, we have two files for our short reads (one for forward reads and one for reverse) and one for long reads. We're going to be referring to these files alot so lets first make a short cut path to these files and then run fastqc on each.

## First, lets check out our sequencing data!

We just got back our high coverage illumina paired-end data and long-read pacbio data from the sequencing facility! Let's see how it turned out. We know that we expect all reads to be 150bp (including adapters) for our short read data, but we're curious _how long our long-read data is_.

First lets specify variables for the paths to our data since we'll be working with them alot
```bash
shortreads="/mnt/data/fastq/shortreads/"
longreads="/mnt/data/fastq/longreads/"
```

Lets break this down into steps.

1) Subset only lines containing read information from our fastqs
```bash
cat $longreads/SalmonSim.Stabilising.p1.3.30k.PacBio.fastq.gz | head #why doesn't this work?
zcat $longreads/SalmonSim.Stabilising.p1.3.30k.PacBio.fastq.gz | grep "chr" -A1 #what does the A option do?
#dont forget about man to read the manual of the command grep!

#we still have some extra information being printed.
#what could you add to this pipe to only keep the sequence? hint: grep -v -E "match1|match2"
```

2) Get the length of every line (read)
```bash
#pipe the output above (every line = seperate read) to the following awk command. 
awk -F "" '{print NR "\t" NF}'  #from the output, can you figure out what NR and NF is? what does the -F "" option do?
#save the output to longread_lengths.txt
```

Instead of just investigating by eye, it would be nice to get some quick stats of the read lengths of these datasets _i.e. what is the average read length for our long read datasets? how much variation is there around the mean?_ \

While R is typically the go to place for statistical analyses, bash can also handle doing some intuitve stats, which will save us the headache of importing data into R.

For example, awk is a super useful tool for working with column and row based analyses. A one-liner in awk can help us calculate the mean read length.

```bash
awk '{ total += $2 } END { print total/NR }' longread_lengths.txt #this gets the mean of column two.
```

**total += $2** <= set the variable "total" equal to the sum of all items in column 2 \

**print total/NR** <= once finished (END), print the column sum after dividing by NR (number of rows) \


We can get the quartiles of read length pretty easily with bash and some simple arithmetic.

```bash
LINECOUNT=`wc -l longread_lengths.txt | cut -d" " -f1`
FIRSTQUART=`echo "$LINECOUNT * 1 / 4" | bc` 
THIRDQUART=`echo "$LINECOUNT * 3 / 4" | bc` 
cut -f2 longread_lengths.txt | sort -n | sed -n "$FIRSTQUART p" 
cut -f2 longread_lengths.txt | sort -n | sed -n "$THIRDQUART p" 
```
**wc -l** <= number of lines \

**cut -d" " -f1** <= keeps only the first column (based on a space delimeter) \

**sed -n "N p"** <= prints (p) the line at value N \


Nice. While theres some variance in our long-read lengths, its nice to see its actually quite consistent. 

Another thing we might want to know is how much coverage our reads give us for each locus in our genome. We are working with a genome size of 10Mb. What is the expected estimate of coverage/bp in just the forward reads?


```bash
READLENGTH=`cat longread_lengths.txt | awk '{ total += $2 } END { print total/NR }' `
NUMREADS=`wc -l longread_lengths.txt | cut -d" " -f1`
MEANCOV=`echo "$READLENGTH * $NUMREADS / 10000000" | bc`
echo "mean coverage = $MEANCOV" 

```

Both the average read length and expected coverage are good things to know for setting expectations for your genome assembly!

## Genome Assembly

Genome assembly can take a long time. Because our course is short, we won't have time to run these commands in tutorial. Instead, today we'll focus on comparing the quality of the assemblies produced through a few different programs, using short, short + long, or long read only data. 


#### Short read assembly: SPADES - *don't run*
#using a single core this takes ~14 minutes
```bash
#mkdir spades
#/ohta1/julia.kreiner/software/SPAdes-3.
#15.3-Linux/bin/spades.py \
#	--pe1-1 ${shortreads}/SalmonSim.Stabilising.p1.1.6400000_R1.fastq.gz \
#	--pe1-2 ${shortreads}/SalmonSim.Stabilising.p1.1.6400000_R2.fastq.gz \
#	-o ./spades/
````


#### Hybrid (short + long) assembly - SPADES & HASLR - *don't run*
#haslr takes ~15 minutes
```bash
#make new dir
#mkdir hybridspades
#run
#/ohta1/julia.kreiner/software/SPAdes-3.15.3-Linux/bin/spades.py \
#	--pe1-1 ${shortreads}/SalmonSim.Stabilising.p1.1.6400000_R1.fastq.gz \
#	--pe1-2 ${shortreads}/SalmonSim.Stabilising.p1.1.6400000_R2.fastq.gz \
#	--pacbio ${longreads}/SalmonSim.Stabilising.p1.3.30k.PacBio.fastq.gz \
#	-t 20 \
#	-o /hybridspades/ 

#mkdir hybridhaslr
#install
#conda install -c bioconda haslr
#run
#haslr.py \
#	-t 20 \
#	-s ${shortreads}/SalmonSim.Stabilising.p1.1.6400000_R1.fastq.gz ${shortreads}/SalmonSim.Stabilising.p1.1.6400000_R2.fastq.gz \
#	-l ${longreads}/SalmonSim.Stabilising.p1.3.30k.PacBio.fastq.gz \
#	-x pacbio \
#	-g 10m \
#	-o hybridhaslr \
#	--cov-lr 40 #takes 15 minutes
#
```

#### Long read assembly - FLYE - *don't run*
#this took 8 minutes with 20 threads

```bash
#install
#conda install flye

#run flye assuming lower quality pacbio reads
#flye \
#	--pacbio-raw ${longreads}/SalmonSim.Stabilising.p1.3.30k.PacBio.fastq.gz \
#	--threads 20 \
#	-o flye/ 
#	--genome-size 10m
```

Now we have four assemblies, 1 short read (spades), 2 hybrid (spades, haslr), and one long-read (flye).

## Assess quality of assemblies

Lets compare assemblies. First, copy these out of the /mnt/data/fasta directory

```bash
mkdir assemblies
cd assemblies

cp /mnt/data/fasta/spades_shortreadonly.fasta ./
cp /mnt/data/fasta/spades_hybrid.fasta ./
cp /mnt/data/fasta/haslr_hybrid.fasta ./
cp /mnt/data/fasta/flye_longread.fasta ./

#we could have done this in one line
#cp /mnt/data/fasta/*.fasta
```

bbmap is command line alignment program that has a collection of nice scripts for library and sequence quality control. We're going to use its stats.sh script to get at some basic stats related to the number and length of sequences in our assembly.

An important stat is N50: the number of contigs making up 50% of your assembly. 
Relatedly, L50 describes for those largest sequences that make up 50% of your assembly, what the minimum sequence length is.
This becomes more intutitive when you look at the lower table outputted by bbmap (unfortunately there is some inconsitency between these definitions in that L/N usage is often swapped, as you'll see later)

We have several assemblies we want to generate stats for, so lets automate this using a for loop:
```bash
mkdir assembly_stats
for i in spades_shortreadonly spades_hybrid haslr_hybrid flye_longread
do

/mnt/software/bbmap/stats.sh in=${i}.fasta > assembly_stats/${i}.stats

done
```

Question: which genome has the best *contiguity* via N50/L50 metrics? Whats an issue with relying just on contiguity?

Having a reference genome of a closely related species can really help asses how *accurate* our assemblies are. In our case we have one better - the high quality reference genome from which our reads were simulated. Lets test out Quast - a program for getting detail stats when comparing to a closely related high quality reference.

Importantly, this allows us to get not just the number and length of our contigs/scaffolds, but also _completeness_ (possible when you know the target genome size) and _correctness_.

```bash
#it would be nice to know how well our assemblies have captured known genes. Lets extract this information from the true reference genome's gff file and format it as quast expects 
head /mnt/data/gff/SalmonAnnotations_forIGV.gff
#quast expects the column order as follows: chromosome, gene id, end, start
awk 'split($9,a,";") {print $1 "\t" a[1] "\t" $5 "\t" $4}' /mnt/data/gff/SalmonAnnotations_forIGV.gff | sed 's/ID=//g' > SalmonReference.genes

#what does $1 and $5 represent in the command above? What about a[1]?

#run quast on all 4 assemblies at once
mkdir quast
/mnt/software/quast-5.0.2/quast.py --help #check out the manual
#notice we're not going to worry about the advanced options 

/mnt/software/quast-5.0.2/quast.py flye_longread.fasta haslr_hybrid.fasta spades_hybrid.fasta spades_shortreadonly.fasta \
	-r /mnt/data/fasta/SalmonReference.fasta \
	-g SalmonReference.genes \
	-o quast

scp -r <username@ip.address>:/home/<usr>/assemblies/quast/ ./
```

Open the report.html results file in your browser and explore the outcome of our assembly efforts. Make sure to eventually click the "View in icarus contig browser".

Question: what new information does this report give us?



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

Printing a particuar row `sed -n "10p" SalmonReference.genes` \

Column means with `awk '{ total += $2 } END { print total/NR }'` \

Assigning variables `shortreads="/mnt/data/shortreads/"` and calling them `echo ${shortreads}` \

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
