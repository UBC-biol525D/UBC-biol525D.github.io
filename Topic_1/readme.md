# Introduction to BIOL525D

## Bioinformatics for evolutionary biologists

*We have tried to keep jargon to a minimum, but if there are things that you want us to clarify or terms you want us to define please don't be afraid to ask!*

In this workshop, we aim to cover the basics of bioinformatics. This course is aimed at those who have some familiarity with computational tools.

There are a number of tools for applying bioinformatic tools with a graphical user interface (such as [Galaxy](https://usegalaxy.org/) and [Geneious](https://www.geneious.com/)). These programs are great and certainly have their place. However, I (Tom) would argue that learning how to use the command line is preferable as it provides far more flexibility and reproducibility and is also a highly transferable skill.

There is a very steep learning curve when it comes to the command line.  The purpose of this workshop is not really to teach you how to use the command line, it takes a lot of practice and learning to get comfortable using the command line. The purpose of this workshop is to demonstrate the fundamentals of bioinformatic analysis. We will need to use the command line throughout the week, so it would be best if participants had taken part in an introduction to the Unix command line workshop before participating in this one, if you not had the time to do that, that's Ok, but there may be places where some things seem a bit opaque. Please do not hesitate to ask questions and do not be put off if you find certain things opaque to begin with, but remember that the primary purpose of this workshop is to become familiar with bioinformatic tools and a general familiarity with UNIX command line environments is strongly encouraged.

Additionally, we have tried to focus on concepts rather than particular software packages. This field is moving so fast that most programs and packages are out of date before too long. That being said, we have had to choose some packages to use for the tutorial, but these should not be seen as the be-all and end-all. There are many packages for specific purposes that we do not have time to go over.

______

## The Integrative Genomics Viewer

In this first tutorial, we are going to use the Integrative Genomics Viewer (IGV) to manually inspect several types of files that you may come across in bioinformatics.

While it may seem topsy-turvy to start the workshop by looking at results before generating results, exploring data in IGV is a really great way to understand the types of data that you will likely encounter in bioinformatics in evolutionary biology.

IGV was developed and is maintained by the Broad Institute (who also maintain many other widely used packages). The IGV provides users with a graphical user interface (GUI) for inspecting and curating datasets. It's a remarkably flexible tool that is invaluable in many instances. Here's how the Broad describes it:

*The Integrative Genomics Viewer (IGV) is a high-performance, easy-to-use, interactive tool for the visual exploration of genomic data. It supports flexible integration of all the common types of genomic data and metadata, investigator-generated or publicly available, loaded from local or cloud sources.*

The first thing that you need to do in this tutorial is to get IGV up and running on your machine.

IGV is written in Java and is available as a pre-compiled package from the Broad Institute. IGV can be freely downloaded at: [https://software.broadinstitute.org/software/igv/](https://software.broadinstitute.org/software/igv/). IGV is written in Java, so if you do not have Java installed on your machine, use the *Java Included* versions of the program for your specific machine.

The download page should look like this:
![](pics/IGV_downloadPage.png)

The links highlighted by the blue blob are what you are after.

## Download data for the tutorial

**Figure out a good way to make the data downloadable - where to host the data so that there is a reproducible ?**

The second thing to do is to download the data package for Tutorial 1 available at the following link:

Once the data has finished downloading, move the data to a memorable location.  

If everything went Ok, you should have the following files:

```shell

ReferenceGenome.fasta.gz ## This is a text file containing the reference genome

SalmonAnnotations.gtf.gz ## A file containing the locations of genomic elements (in this case genes)
SalmonAnnotations.gtf.gz.tbi ## An index for the above file  

Salmon.HiSeq.30x.bam ## A file containing the alignments of paired-end
                     ## Illumina HiSeq reads to the reference genome - at 30x

Salmon.HiSeq.10x.bam ## A file containing the alignments of paired-end
                     ## Illumina HiSeq reads to the reference genome - at 10x

Salmon.ddRAD.bam ## A file containing the alignments of double-digest RAD seq.
                 ## reads to the reference genome

Salmon.Nova.bam ## A file containing the alignments of PacBio NovaSeq reads to
                 ## reads to the reference genome

Salmon.MiSeq.10x.vcf # A file containing variants called from the 10x Illumina data

Salmon.MiSeq.30x.vcf # A file containing variants called from the 30x Illumina data

## Each of the files ending in ".bam" also have an index file (those ending in ".bai")

SalmonData_checksums.sha ## See below

```
## Check data integrity

The following is optional, only do this if you are comfortable working on the command line already - return to this at a later stage if not.

<details>
  <summary>Click to expand!</summary>

Often when downloading large files, we want to perform a check on the integrity of the data downloaded. Perhaps you were downloading a large file and the power went down before you were not able to check the data. More generally, when downloading many large files you may not have the time (or the patience) to manually curate each one to ensure that it was downloaded correctly.

It is good practice to check data integrity when moving files from place to place and there are useful functions for checking data integrity. The two main methods that are used are ```shasum``` and ```md5```. ```shasum``` is perhaps a bit more common as it is available as standard on MacOS. Both methods generate what is called a "checksum", a hexadecimal string (e.g. "d241941bac307bd853fd21945d029e62c83cea71") that is unique to a given file.

If you inspect the contents of ```SalmonData_checksums.sha``` you'll notice that there are two columns in the file, the left hand column contains a bunch of checksums, the righthand column contains the name of corresponding files. You can compare all the files in the directory you downloaded using:


**Store the checksums with the data, whereever we host it.**

```shell
shasum -c SalmonData_checksums.sha

```
</details>


## Load the reference genome into IGV
*Let us know if you've had difficulties setting up IGV on your machine*

Open up IGV on your machine, you should be seeing something like the following (don't worry if it looks a little different):
![](pics/IGV_startScreen.png)

You may notice "hg19" in the drop down menu on the top left (highlighted in orange). "hg19" stands for human genome version 19, it's fun to explore the human genome, but today we're going to explore the genome of the system we're working with in this workshop .

The first thing to do is to load in our reference genome. From the drop-down menu in the top left, choose the ```Genomes``` drop-down menu and choose ```Load Genome From File...```. When the box opens up, navigate to the data you downloaded and choose the ```ReferenceGenome.fasta.gz``` file. The file extension ".fasta.gz" tells us that this is a FASTA file (a simple text file that represents a genetic sequence - DNA, RNA or peptides).

Once you've done that, you should be able to select "SalmonReference.fasta" from the drop down menu. Explore IGV and figure out how many chromosomes our reference genome contains and how long the chromosomes are.



# Load sequence annotations into IGV



# Load sequence alignments into IGV

### Low coverage Illumina sequencing

Let's go ahead and load a file containing sequence alignments. To begin with, load up the file ```Salmon.MiSeq.10x.bam``` in IGV. You do this by selecting "SalmonReference.fasta" from the dropdown menu (where you can see)

This file (and the others like it) contains information on paired-end short reads generated using an Illumina MiSeq machine. The data has been aligned to the reference genome we are using and is stored as a "BAM" file. A "BAM" file is a file containing the alignment of genetic sequences to a reference genome in binary format. Files in binary format are often smaller than raw text, and going from binary to text involves work for your computer, so binary can be worked with more efficiently. BAM stands for **B**inary **A**lignment **M**ap. Don't worry, we'll go over the details of SAM/BAM formats later in the week, for now though, just know that they are a type of file containing genetic sequences.

If the data loaded happily, choose a chromosome from the dropdown menu and zoom in to any location in the genome (by double clicking or using the zoom bar in the top right corner). The window should look a little like this:
<details>
  <summary>Click to expand!</summary>
  ![](pics/IGV_HiSeq_lowCoverage.png)
</details>


We are looking at a graphical representation of the short reads aligned to our refernece genome. Illumina technology is amazing, but it is not perfect. Explore the alignment, using your mouse or trackpad, scroll through the alignment and get a feel for navigation in IGV.

#### Discussion points:

* *What do you think the challenges are with using this data to make biological inferences?*
* *What features should be used to assess the quality of an alignment?*




### Higher coverage Illumina sequencing

After you've had a little play with the reference genome, let's go ahead and load a file containing sequence alignments. To begin with, load up the file ```Salmon.MiSeq.10x.bam``` in IGV. You do this by selecting "SalmonReference.fasta" from the dropdown menu (where you can see)

This file (and the others like it) contains information on paired-end short reads generated using an Illumina MiSeq machine. The data has been aligned to the reference genome we are using and is stored as a "BAM" file. A "BAM" file is a file containing the alignment of genetic sequences to a reference genome in binary format. Files in binary format are often smaller than raw text, and going from binary to text involves work for your computer, so binary can be worked with more efficiently. BAM stands for **B**inary **A**lignment **M**ap. Don't worry, we'll go over the details of SAM/BAM formats later in the week, for now though, just know that they are a type of file containing genetic sequences.

If the data loaded happily, choose a chromosome from the dropdown menu and zoom in to any location in the genome (by double clicking or using the zoom bar in the top right corner). The window should look a little like this:

![](pics/IGV_startScreen.png)


## I hope you had fun!
