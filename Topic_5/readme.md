---
title: "Topic 5: Sequence Alignment"
permalink: /Topic_5/
topickey: 5
topictitle: "Sequence Alignment"
---

### Accompanying material

* [Slides](./Topic 4.pdf)


Today we're going to align sequence data to a reference genome using BWA and explore what a BAM file is.

Let's set up a directory structure so the resulting files will be organized and copy the raw data to your home directory.

```bash

# Navigate to your working directory
cd ~

# Copy the reference genome to your working directory
cp -r /mnt/data/fasta ./

# Copy the fastq files to your working directory
cp -r /mnt/data/fastq/GWAS_samples/ ./

# Make a new directory for your resulting BAM files
mkdir bam
```
We are going to work with the true genome of the species that you explored yesterday because we have limited time. However, can you think of how the choice of assembly would affect the mapping of our data?

First let's index our reference genome.

```bash
# Index the references for BWA.

bwa index fasta/SalmonReference.fasta

```

Now finally we can run BWA and align our short read data.

Using the approach 1 genome assembly:


```bash

bwa mem \ # We have installed BWA on the VMs, but that might not always be the case
  fasta/SalmonReference.fasta \
  GWAS_samples/Salmon.p1.3.i1.400000_R1.fastq.gz \
  GWAS_samples/Salmon.p1.3.i1.400000_R2.fastq.gz \
  -t 2 \
  -R '@RG\tID:sample_1\tSM:1\tPL:illumina\tPU:biol525d\tLB:sample_1_library' \
  > bam/Salmon.p1.3.i1.sam

```
*This will take a few moments to run*


Lets break this command down since it has several parts:
**/usr/bin/bwa** <= We're calling the program _bwa_ from the directory _/usr/bin/_. This is the full path to that program so you can call this no matter where you are in the file system.

* **mem** <= This is the mode of bwa we are running. It is an option specific to bwa and not a Unix command.

* **\\** <= Having this at the end of the line tells the shell that the line isn't finished and keeps going. You don't need to use this when typing commands in, but it helps break up really long commands and keeps your code more organized.

* **fasta/SalmonReference.fasta** <= This is the reference genome. We're using a relative path here so you need be in /mnt/<USERNAME> or it won't be able to find this file.

* **GWAS_samples/Salmon.p1.3.i1.400000_R1.fastq.gz** <= This is the forward read (e.g. read 1)  set for the first sample. It's also a relative path and we can see that the file has been gzipped (since it has a .gz ending).

* **GWAS_samples/Salmon.p1.3.i1.400000_R2.fastq.gz** <= This is the reverse read (e.g. read 2)  set for the first sample.

* **-t 2** <= This is telling the program how many threads (i.e. cpus) to use. In this case we're only using two because we're sharing the machine with the other students.

* **-R '@RG\tID:sample_1\tSM:1\tPL:illumina\tPU:biol525d\tLB:sample_1_library'** <= This is adding read group information to the resulting SAM file. Read group information lets other programs know what the sample name along with other factors. It is necessary for GATK to run later on.

* **> bam/Salmon.p1.3.i1.sam** <= This is directing the output of the program into the file bam/Salmon.p1.3.i1.sam

We now have our reads aligned to the genome in a human readable format (SAM) instead of binary format (bam) which we will use later. Generally we keep our data in BAM format because its more compressed but we can use this opportunity to better understand the format.


Lets examine the SAM file. It contains all the information on the reads from the fastq file, but also alignment information.

```bash

# Let's view that SAM file
less -S bam/Salmon.p1.3.i1.sam

# Notice the @PG line that includes the program call that created the SAM file.
# This is useful for record keeping.

```
### *Note*
The option `-S` when running less chops lines that are longer than the page. This is normally just an aesthetic choice. When looking at SAM/BAM files this is quite necessary!

 

### Questions:
1. How are reads ordered in the SAM file?
2. What does the 6th column represent? What would the string "3M1I3M1D5M" tell you?
3. What are three possible reasons why mapping quality could be low for a particular read?

____________________________

At this point we'll introduce a very useful - and incredibly widely used - piece of software called `samtools`. As the name suggests, `samtools` is a program for working with SAM/BAM files.

### *Note*
`samtools` can produce very useful summaries of alignments - try running `samtools flagstat bam/Salmon.p1.3.i1.sam`.

A question that you might ask of an alignment would be, what proportion of my reads mapped to the genome? At this stage, our SAM file contains all the read data, whether reads mapped or not. Using `samtools`, we can easily get a count of the number of reads that successfully mapped to the genome.


```bash

samtools view -c bam/Salmon.p1.3.i1.sam

```

Lets break this command down:
* `samtools`  - the program tat we want to run
* `view` - the mode we want to run the program in
* `-c` - this flag indicates that we want a count of reads
* `bam/Salmon.p1.3.i1.sam` - The input file


This should have printed the total number of mapped reads to screen. There should be no surprises here.  

Now what we're going to do is to remove the `-c` option, which causes `samtools` to send the output straight to STDOUT. This is handy as it means we can pipe it into another process.

In the following, we'll take our SAM file (human readable) and convert to a BAM file (machine readable) and sort reads by their aligned position.

```bash
samtools view -bh bam/Salmon.p1.3.i1.sam | samtools sort > bam/Salmon.p1.3.i1.sort.bam
```


Lets break this command down:
* `samtools`  - the program tat we want to run
* `view` - the mode we want to run the program in
* `-bh` - this is actually two flags, one that tells samtools to express the data in binary form and the other that tells samtools to include the header
* `bam/Salmon.p1.3.i1.sam` - The input file
* `|` - The pipe symbol - you should be familiar with this by now
* `samtools sort` - another mode of samtools that sorts SAM/BAM files by coordinate
* `> bam/Salmon.p1.3.i1.sort.bam` - the name you want to give to the output file


With this command we're using the pipe "|" to pass data directly between commands without saving the intermediates. This makes the command faster since its not saving the intermediate file to hard disk (which is slower). It can be more risky though because if any steps fails you have to start from the beginning.

Next we want to take a look at our aligned reads. First we index the file, then we use samtools tview.

```bash
samtools index bam/Salmon.p1.3.i1.sort.bam # build an index of a BAM file
samtools tview bam/Salmon.p1.3.i1.sort.bam  --reference ref/approach_1.fasta

#use ? to open the help menu. Scroll left and right with H and L.
#Try to find positions where the sample doesn't have the reference allele.

```

`samtools tview` is similar to IGV but is accessible directly from the command line.

You can jump to a specific location in a BAM file with `samtools tview` using the following command:

```

samtools tview bam/Salmon.p1.3.i1.sort.bam  --reference fasta/SalmonReference.fasta -p chr_1:80000 

```
The additional option ` -p chr_1:80000 ` tells `tview` to jump straight to chr_1 position 80000. Any valid location in the SAM/BAM can be referenced that way.

### *Note*

Another useful summary that `samtools` can produce very quickly is coverage stats. Try running `samtools depth bam/Salmon.p1.3.i1.sort.bam`. Can you think of how you could use the tools you were learning yesterday could be used to take the output from `samtools` to quickly calculate the average depth?



## Exercise


A bash script is a plain text file (i.e. not rich text, nor a word doc) which contains bash commands. You can create the file on your computer and copy it over, or you can edit it directly from the server with one of the installed editors (this is covered in [topic 2, Editing](../Topic_2/#editing). The name of the file is up to you, but bash scripts are given the `.sh` extension by convention.

Here's an example of what you might see inside a bash script:

```
age=10
echo "I am $age years young"

```
(we used a similar example the other day)

The lines of this script do the following:

* `age=10` - this assigns the number 10 to a variable called `age`
* `echo "I am $age years young"` - this prints a piece of text to the screen containing that age variable



So a bash script (a.k.a. shell script) is a just set of commands saved as a file. If you named the shell script from above as a file named `myScript.sh`, you could execute the commands by simply running:

```bash
sh myScript.sh
```

_____________________________

Obviously that example is a little silly, but hopefully you can see how writing shell scripts is a very useful and efficient way of organising your work at the command line. 

If you look in the `/mnt/data/fastq/GWAS_samples/` directory, you'll see that we have data here for 10 samples. It would be very tedious to align each one of these as we have for the single file above.

For this exercise, try writing a bash script to produced a sorted BAM file for each sample as you have done for the single sample above. 

When writing a shell script, try to think of the steps that do not need to be repeated over and over again.


HINTS:
  * Use variables for directory paths e.g. `bwa=/mnt/bin/bwa-0.7.17/bwa`
  * Use a loop.
  {: .spoiler}

MORE HINTS:
  * for/while loops can receive input from stdin (or a file):

        while read fname;
	  do echo processing "$fname";
	done < list_of_things.txt
  {: .spoiler}

  * You can do pathname manipulation with `basename` and `dirname` (see manual pages):

    ```
    dirname a/b/c          # prints a/b
    basename a/b/c.gz      # prints c.gz
    basename a/b/c.gz .gz  # prints c

    fpath=/path/to/the/file.gz
    base=$(basename "$fpath")    # assign "file.gz" to variable base
    echo "$base"                 # prints file.gz
    ```
  {: .spoiler}

<details>
<summary markdown="span">**Answer**
</summary>
```bash
  # First set up variable names
  # These may be slightly different on the VMs
  bam=~/bam
  fastq=~/GWAS_samples
  bwa=bwa
  ref_file=~/fasta/SalmonReference.fasta

  #Then get a list of sample names, without suffixes
  ls $fastq | grep R1.fastq.gz | sed s/_R1.fastq.gz//g > $bam/samplelist.txt

  #Then loop through the samples
  while read name
  do
    $bwa mem \
    -R "@RG\tID:$name\tSM:$name\tPL:ILLUMINA" \
    $ref_file \
    $fastq/${name}_R1.fastq.gz \
    $fastq/${name}_R2.fastq.gz \
    -t 1 > $bam/$name.sam;

    samtools view -bh $bam/$name.sam |\
    samtools sort > $bam/$name.sort.bam; 
    samtools index $bam/$name.sort.bam
    
    rm $bam/$name.sam # Remove intermediate file
    
  done < $bam/samplelist.txt
```
</details>

After your final BAM files are created, and you've checked that they look good, you should remove intermediate files to save space. You can build file removal into your bash scripts as we've done in the worked example, but it is often helpful to only add that in once the script is up and running. It's hard to troubleshoot a failed script if it deletes everything as it goes.

### By topic 7, you should have created cleaned BAM files for all samples.

## Questions for Discussion
1. Is an alignment with a higher percent of mapped reads always better than one with a lower percent? Why or why not?
2. I want to reduce the percent of incorrectly mapped reads when using BWA. What setting or settings should I change in BWA?
3. What are two ways that could be used to evaluate which aligner is best?
