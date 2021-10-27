---
title: "Topic 6: RNAseq analysis exercise"
permalink: /Topic_6/
topickey: 6
topictitle: "RNASeq Analysis"
---

In this tutorial we are going to do the following:

1. Align the RNAseq data to the reference
2. Obtain counts of the reads mapping to genes in the Salmon genome
3. Test for differential gene expression between fish from warm and cold parts of the river - This last step will be done in Rstudio

We have data for 12 fish, 6 from each location in the river.

FASTQ files containing RNA seq reads for the 12 fish are located at:

```
/mnt/data/fastq/rna/
```

We've got paired end reads, so there's two fastq files per sample. The samples are named in a pretty obvious way, but I'll let you see that for yourself.

# 1. Aligning RNAseq data

We will use the splice-aware RNA read aligner STAR. As mentioned in the lecture, STAR performs well under default parameter settings. For this tutorial, we will not have time to tinker with alignment settings, we will just use the program defaults.

## Install and locate relevant software

We have preinstalled STAR on each of the servers, but if you were to do this yourself you can follow the easy instructions on the STAR GitHub page ([https://github.com/alexdobin/STAR](https://github.com/alexdobin/STAR)).

STAR is in the following location on each of the VMs:
```sh
/mnt/software/STAR-2.7.9a/source/STAR
```

It is a bit annoying to have to type that whole path in all the time, so let's add it to our PATH.

The PATH is a list of directories that your shell searches through when you enter a command.

If you add the directory containing the STAR executable to the path, you'll be able to access it at the command line without the need to specify the full path.

Here's how you can do that:
```sh
# navigate to home directory
cd

# open up your bash profile
emacs .bashrc  # I like emacs, but use whichever text editor you prefer

```

When in your text editor, add the following text and save the file:
```
echo "Hello, you are cool!"
export PATH="/mnt/software/STAR-2.7.9a/source/:$PATH"
```

To get your operating system to work with the new instructions you've given it you can refresh your terminal session using the following command:

```
source ~/.bashrc
```
This command just reloads the configuration file you've just edited. You can also accomplish this by logging out and logging back in.

There are a couple of things here that we have not covered so far. One is the shell profile. This is a set of commands that the server runs when you log in. As you develop your skills it becomes really handy to customise your bash profile.

Now check that it worked:

```
STAR --help

```
That should access the STAR executable and print a bunch of help text to the screen.

## Build a genome index

Now that we have STAR up and running, the first thing we'll need to do is build an index for the reference genome so that we can align our reads to it. Building an index is necessary for many bioinformatic operations. It is not dissimilar to the index of a book that tells you where to find certain topics or key wrods. In the case of a genomic index, one builds a record of where you find certain sequences. That makes looking through the genome much more efficient.


STAR, like many read aligners, has many modes of operation and building an index is just one of them. To tell STAR

```sh
STAR --runThreadN 2 \ # The number of threads to spawn this process on
                              --runMode genomeGenerate \ # The mode of operation for STAR
                              --genomeDir location_to_save_index \ # A place to store the index file
                              --genomeFastaFiles SalmonReference.fasta \ # The location of the reference genome
                              --sjdbGTFfile SalmonAnnotations_forIGV.gff # The location of the genome annotations, in GFF format


```

You'll need to specify the location of the Salmon reference genome and genome annotation files appropriately.

This takes about half a minute to run on 2 threads. Each VM has a maximum of 16 threads, so please do not use too many at once!


## Map RNA-seq reads

Now we've built our index, we can map our reads to the reference genome.

Let's use the STAR to actually map the reads:

```sh

STAR --genomeDir location_to_save_index/ \ # This tells STAR where we've put the reference genome - the place you specified above
      --readFilesIn cold_sample_01_1.fq.gz cold_sample_01_2.fq.gz \ # Give the two FASTQ files for paired-end reads
      --outFileNamePrefix cold_sample_04. \ # Give a prefix for all of the output files
      --outSAMtype BAM SortedByCoordinate \ # This tells STAR to outut the alignments in BAM format and sorted by coordinate
      --outSAMunmapped Within \ # Puts the unmapped reads into the BAM file
      --outSAMattributes Standard \ # Use standard SAM formatting
      --readFilesCommand zcat \ # This tells STAR that the fastq files were gzipped
      --runThreadN 2

```
You can use multiple threads when aligning reads with STAR too. When you're working with full sized datasets that will for sure come in handy, but with the data we are working with that's not a big issue. Mapping the reads using a two threads takes about a minute per sample .

### Challenge 1
The above mapped reads for a single sample. We need to repeat the above but for all samples that we have data for. Can you think of a way to run the above for all the samples we have data for that does not involve manually typing each one in?

Once we have BAM files for each sample we'll need to count the number of reads associated with each gene.

# Obtain raw read counts

To obtain counts of reads we will use the `htseq-count` tool. This is a popular and very handy tool that can be used to obtain counts of RNA seq reads that have been mapped to a genome.

`htseq-count` sends its list of read counts to STDOUT, so you need to capture it with a redirect if you want to save the output
```
htseq-count -s no \
            -r pos \
            -t exon \ # What type of feature will our data have mapped to?
            -i gene \
            -f BAM \
            cold_sample_07.Aligned.sortedByCoord.out.bam \
            /mnt/data/anno/SalmonAnnotations.gff  > cold_sample_07.read_counts.txt
```

The file
