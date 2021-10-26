---
title: "Topic 3: Bioinformatics Gotchyas"
permalink: /Topic_3/
topickey: 3
topictitle: "Bioinformatic Gotchyas"
---

## Accompanying material
[Lecture Slides](./Topic_3.pdf)


# Checking data integrity

Often when downloading large files, we want to perform a check on the integrity of the data downloaded. Perhaps you were downloading a large file and the power went down before you were not able to check the data. More generally, when downloading many large files you may not have the time (or the patience) to manually curate each one to ensure that it was downloaded correctly.

It is good practice to check data integrity when moving files from place to place and there are useful functions for checking data integrity. The two main methods that are used are shasum and md5. shasum is perhaps a bit more common as it is available as standard on MacOS. Both methods generate what is called a “checksum”, a hexadecimal string (e.g. “d241941bac307bd853fd21945d029e62c83cea71”) that is unique to a given file.

The data that you obtained for Topic 1 also contained a file called ```SalmonData_checksums.sha```. If you inspect the contents of ```SalmonData_checksums.sha``` you’ll notice that there are two columns in the file, the left hand column contains a bunch of checksums, the righthand column contains the name of corresponding files. You can compare all the files in the directory you downloaded using:

Store your checksums with your data, whereever you host it.

Navigate to the location of the data and check that the downloaded data actually match those that were sent by the  
```
shasum -c SalmonData_checksums.sha
```
# Inspect read quality

## Lets look at some properties of our reads
```bash
shortreads="/mnt/data/fastq/shortreads/"
longreads="/mnt/data/fastq/longreads/"
mkdir assembly
cd assembly
#fastqc is a widely used program for looking at read quality
fastqc ${shortreads}/*fastq.gz ${longreads}/*fastq.gz -o ./ 
#some programs allow you to input multiple files at once using wild cards. this is convienient. 
#otherwise we would have had to specify each file individually. heres anotherway we could use wildcards here: 
#fastqc ${shortreads}/*R?.fastq.gz ${longreads}/*fastq.gz -o ./ )
#summarise into a single report
multiqc ./
```
Download the output of multiqc run as follows:
```bash
#From a terminal window running on your local computer
scp <username@ip.address>:/home/<username>/assembly/multiqc_report.html <path on your computer where you want the file>
```
Since we are using simulated data, our data does not have adapters. But taking a look at the signatures of our highly quality reads can help you get a sense of what very high quality data might look like, and additionally,by comparing to your real data, how sequencing technology influences data quality.  
Discussion Quesiton: our data effectively has had adapters removed from reads and already been trimmed for low quality BPs near the end of reads. what might be the cons of read trimming? 



# Exercise: Convert GFF to BED

A routine task in bioinformatics is to convert files from one form to another. In this exercise we'll use the Unix command line to take a file from one format and make the necessary changes to put it in another.

For this exercise, pretend that a collaborator of yours has written a piece of analysis software. This analysis examines the total length of the coding portion of the genome and does a sophisticated evolutionary analysis. However, your colleague was not very helpful and did not know that the data is currently in GFF format. The program they wrote only takes BED files. Here's the specification of the file format from your collaborator:

_ _ _ _ _ _ _ _ _ _
*My new and improved algorithm for analysing selection in the genome requires, as input, a BED file in the following format:*
* *One file per chromosome*
* *Chromosomes are referred to using a single integer name*
* *Only the coordinates of the coding sequences (CDS)*
_ _ _ _ _ _ _ _ _ _

For this exercise can you 

