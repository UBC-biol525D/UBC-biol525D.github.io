---
title: "Topic 3: Bioinformatics Gotchyas"
permalink: /Topic_3/
topickey: 3
topictitle: "Bioinformatic Gotchyas"
---

## Accompanying material
[Lecture Notes](https://drive.google.com/file/d/1EkfhMAsj9vgr2wmgD1qrTOKgE7TKFArX/view?usp=share_link)

# Part 1: Overview of sequence data

In this exercise, let's pretend that you have been notified that your data is ready for download from the sequence centre. You submitted a short read sequencing library and a long read sequencing library generated using Illumina and PacBio technologies, respectively.

To start, make a local copy of the data files in your home directory. Run the following:
```
cp /mnt/data/Tutorial_3_data.tar.gz ./ # Make a local copy of the sequence data

tar -zxvf Tutorial_3_data.tar.gz  # Extract the data from the tarball that it was shipped in

```

A "tarball" is a way of compressing an entire directory full of files so that you only have to download one data packet rather than many. The options given to tar do the following:
* ``-z``, this tells the program that the file "Tutorial_3_data.tar.gz" is compressed or zipped
* ``-x``, this tells the program to extract the files from the archive
* ``-v``, this tells the program to print a log to screen
* ``f``, this tells the program that you are specifying the tar-ball as a file on the command line


## Checking data integrity

Often when downloading large files, we want to perform a check on the integrity of the data downloaded. Perhaps you were downloading a large file and the power went down before you were not able to check the data. More generally, when downloading many large files you may not have the time (or the patience) to manually curate each one to ensure that it was downloaded correctly.

It is good practice to check data integrity when moving files from place to place and there are useful functions for checking data integrity. The two main methods that are used are shasum and md5. shasum is perhaps a bit more common as it is available as standard on MacOS. Both methods generate what is called a “checksum”, a hexadecimal string (e.g. “d241941bac307bd853fd21945d029e62c83cea71”) that is unique to a given file.

The data that you obtained for Topic 1 also contained a file called ```SalmonData_checksums.sha```. If you inspect the contents of ```SalmonData_checksums.sha``` you’ll notice that there are two columns in the file, the left hand column contains a bunch of checksums, the righthand column contains the name of corresponding files. You can compare all the files in the directory you downloaded using:

Store your checksums with your data, where ever you host it.

Navigate to the location of the data and check that the downloaded data actually match those that were sent by the  
```
cd Tutorial_3_data
shasum -c SalmonData_checksums.sha
```

As long as everything worked well, you should have seen reassuring messages print to screen after invoking the `shasum` command.

______________________________

Now we've got our data, let's just make a couple of directories to tidy things up.

```
mkdir short_reads
mv Salmon.Illumina.R?.fastq.gz short_reads/ # Note the use of the single character wildcard!

mkdir long_reads
mv Salmon.PacBio.fastq.gz long_reads/

cd ../

```

Great, that's nice and tidy now.

## How many reads do we have?

A question that might be on top of mind would be, how many reads are we working with?

Can you use command line tools to get a count of the number of reads that we have from the two datasets?

*Hint, remember the structure of a fastq file*


## Inspect read quality

### Install fastqc

We are going to use a program called `FastQC` to examine the quality of the sequencing reads that we are working with.

It is really easy to generate a report using `FastQC`, it's just a matter of telling the program where the files are that you want to analyse. For example:

```bash

shortreads="/home/tommyB/Tutorial_3_data/short_reads/" # Obviously use your own user name here
longreads="/home/tommyB/Tutorial_3_data/long_reads/"

mkdir qualityControl
cd qualityControl

```

Most programs won't be install for you when you start analyzing your data. Its good practice to get a sense of how to do this!

The first thing you need is the address of the most recent version of `fastqc`. If you navigate the program's website, you'll find several download links when you click "Download Now":
[https://www.bioinformatics.babraham.ac.uk/projects/fastqc/](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)

*DO NOT LEFT CLICK ON THE DOWNLOAD LINK!*
Instead, right click the link for Linux downloads and hit "Copy link address"

Now, we can use that link to download the package directory to our local directories:
```bash
wget <insert copied link here>

# This should have downloaded a zipped file called something like fastqc_v0.11.9.zip

# To unzip just run:
unzip fastqc_v0.12.1.zip

# This will have made a directory, go into it:
cd FastQC

## Make the fastqc file executable:
chmod 755 fastqc ## Add execute permissions for all users

```

Phew, now we should have a happy version of fastqc that we can work with.

### Now lets look at some properties of our reads. We're including read sets generated from different sequencing platforms resulting in different read lengths (Illumina/Pacbio), but also short reads from historical specimens. We need to make a path to the historical read set, like we did above.

### We can also add a path to programs, so we can just type the program name

```bash

historical_shortreads=~/Tutorial_3_data/historical_shortreads/
fastqc=/home/ubuntu/Tutorial_3_data/qualityControl/FastQC/fastqc

#OK lets run the quality control analysis now that we're all set
cd ~/Tutorial_3_data/qualityControl

#fastqc is a widely used program for looking at read quality
#You might be wondering how to run fastqc? How would you figure it out on your own? Remember that pretty much all programs have a manual, which can often be called by appending --help to the program name. Usually the most important information is at the top i.e. the recipe for running the program.

$fastqc --help

#You'll note that fastqc is conveniently coded to allow inputting multiple input files at once, this isn't the case for alot of programs. Make sure to get into the habit of checking out the usage section of a program's manual before you run it.

#OK lets run fastqc
$fastqc  ${shortreads}/*fastq.gz ${historical_shortreads}/*fastq.gz ${longreads}/*fastq.gz -o ./

#Here's another couple of ways we could use wildcards here
##$fastqc ${shortreads}/*R?.fastq.gz ${historical_shortreads}/*R?.fastq.gz ${longreads}/*fastq.gz -o ./ )
##$fastqc ${shortreads}/*R[1,2].fastq.gz ${historical_shortreads}/*R[1,2].fastq.gz ${longreads}/*fastq.gz -o ./ )

#summarise into a single report
multiqc ./
```
Download the output of multiqc run as follows:

```bash
#From a terminal window running on your local computer
scp <username@ip.address>:/home/<username>/assembly/multiqc_report.html <path on your computer where you want the file>
```

Open that file up in an internet browser and have a look.


**Discussion: The "shortread" and "longread" data is simulated (with low error rates), whereas the "historical_shortread" is actually from DNA prepared from 100 year old samples. Historical sample preparation tends to use PCR amplification to maximize DNA yield and might show signals of degradation. Do you see signals of this?**



# Part 2: Convert a GFF file to BED format

A routine task in bioinformatics is to convert files from one form to another. In this exercise we'll use the Unix command line to take a file from one format and make the necessary changes to put it in another.

For this exercise, pretend that a collaborator of yours has written a piece of analysis software. This analysis examines the total length of protiein-coding genes and does a sophisticated evolutionary analysis with that information. However, your colleague was not very helpful and did not know that the data is currently in GFF format. The program they wrote only takes BED files. Here's the specification of the file format from your collaborator:

_ _ _ _ _ _ _ _ _ _
*My new and improved algorithm for analysing selection in the genome requires, as input, a files with the following characteristics:*
* *One file per chromosome*
* *Chromosomes are referred to using a single integer name*
* *The input **MUST** be in BED format*
* *The program can detect when there is duplicate sequences, but it is **much** faster if it does not have to identify them*
_ _ _ _ _ _ _ _ _ _

For this exercise can you generate BED files that your collaborator's software would accept?

The annotation file for the entire genome is on the server - we suggest you make a local copy and store it in a convenient location. Here's a command to grab the annotations:

```
cp /mnt/data/anno/SalmonAnnotations.gff ./
```

### If you don't know the commands, call one of the instructors over to walk through the steps you would take to complete the exercise



<details>
<summary markdown="span">**Click for partial solution!**
</summary>
```bash
   > cat /mnt/data/anno/SalmonAnnotations.gff | awk 'BEGIN {OFS = "\t"};{if ($3=="gene") print  $1,$4-1,$5}'

#this gets you part of the way there, but the chromosome values (column 1) still aren't single integer. the command `tr` (stands for translate) could be useful for this, or `sed 's/find/replace/g'`

# The program `bedtools` would come in handy to deal with duplicate entries

```
</details>
