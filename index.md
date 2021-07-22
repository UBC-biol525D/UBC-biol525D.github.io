---
title: "Biol525D Homepage"
layout: home
menuItem: "General info"
menuPosition: 1
---

<h1>{{ site.courseName }}</h1>

<img src="{{ site.baseurl }}/style/header.jpg" width="100%" alt="Nitobe Memorial Gardens -- CC BY-NC 2.0 - @kardboard604 https://www.flickr.com/photos/moov4/10157685024/" title="Nitobe Memorial Gardens -- flickr @kardboard604">

## Description
The purpose of this course is to provide graduate students with the theoretical knowledge and practical skills for the evolutionary analysis of next generation sequence data. The course will entail data retrieval and assembly, alignment techniques, variant calling, gene expression analyses, hypothesis testing, and population genomic and phylogenomic approaches. The course will be presented as a series of short lectures and lab exercises over a one week period in August.

## Instructors
Dr. Kathryn Hodgins, Dr. Tom Booker, Julia Kreiner

## Format
A mix of lecture and lab exercises, running in a 2-hour block.

## Prequisites
Students should have basic knowledge in R and some command line knowledge (although the latter could be obtained during the course)

## Evaluation
Participation in discussions and lab exercises.

Consult [Grading Information](resources/grading.md).

## Assignments
Consult [Daily Assignments](resources/daily_assignments.md)

## Lectures/seminars
30th November 2020 - 4th of December 2020.

There are two sessions a day, from 10am to Noon and from 1pm to 3pm.

The course will be run entirely over Zoom - links to the sessions will be sent to those who registered interest. 

## Basic course structure

The course material is organized in several topics, with slides and coding examples.

To get up to speed on working with a Unix system, take a look at the [unix help](resources/unix_ref.pdf) file. There are some resources there that will help you find the specific command you need for each task.

## Syllabus
1. [Topic 11](./Topic_11/) Broad introduction: Scope of course, goals, overview of technology and bioinformatics, and the future of sequencing [Booker,T]
2. [Topic 2](./Topic_2/) Programming for biologists [Booker, T]
3. [Topic 3](./Topic_3/) Fastq files and quality checking/trimming [Kay, H]
4. [Topic 4](./Topic_4/) Alignment: algorithms and tools [Booker, T]
5. [Topic 5](./Topic_5/) Assembly: transcriptome and genome assembly [Kay, H]
6. [Topic 6](./Topic_6/) RNAseq + differential expression analysis [Kay, H]
7. [Topic 7](./Topic_7/) SNP and variant calling [Kreiner, J]
8. [Topic 8](./Topic_8-9/) Population genomics and plotting in R (Part 1) [Kreiner, J]
9. [Topic 9](./Topic_8-9/) Population genomics and plotting in R (Part 2) [Kreiner, J]
10. Topic 10 Case studies [Kreiner,J/Booker,T]

In previous years, there has been a short introduction to phylogenetic analysis included in the workshop. We have not included that topic this year, but the slides and materials are available for those who are interested:

11. [Topic 11](./Topic_10/) Phylogenetic analysis


## Obtaining all the files on this site

You may use your internet connection to browse this site, or
you may download the entirety of the files on the site in one
constantly updated zip archive
[here](https://github.com/owensgl/biol525D/archive/master.zip)

This method dosesn't require `git`, however, you'll have to manually
update the files this way (by downloading the whole repo again).

To obtain to all the files via git, type:

    git clone https://github.com/UBC-biol525D/UBC-biol525D.github.io.git

To update the all the files at any point in the future, change to the **biol525D** directory that was created by the previous command and type:

    git pull


## Use and modification of these resources

You may use any of the materials provided here, and modify them in any way, provided there is appropriate attribution according the license found below and included with this project.

## License and Copyright

Copyright (C) 2015 S. Evan Staton, Sariel Hubner, Sam Yeaman

Modified work (c) 2016, 2017, 2018 Gregory Owens, Kathryn Hodgins

Modified work (c) 2019 Gregory Owens, Kathryn Hodgins, J.S. Legare

Modified work (c) 2020 Kathryn Hodgins, Julia Kreiner, Tom Booker

Modified work (c) 2021 Tom Booker, Julia Kreiner

This program is distributed under the MIT (X11) License, which should be distributed with the package.
If not, it can be found here: http://www.opensource.org/licenses/mit-license.php

## About this site:

   This site is powered by GithubPages, and the code backing it is on GitHub [here](https://github.com/UBC-biol525D/UBC-biol525D.github.io).
