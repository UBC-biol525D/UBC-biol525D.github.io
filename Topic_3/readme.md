---
title: "Topic 3: Bioinformatics Gotchyas"
permalink: /Topic_3/
topickey: 3
topictitle: "Bioinformatics Gotchyas"
---

## Accompanying material
[Lecture Slides](./Topic_3.pdf)

## Checking data integrity

Often when downloading large files, we want to perform a check on the integrity of the data downloaded. Perhaps you were downloading a large file and the power went down before you were not able to check the data. More generally, when downloading many large files you may not have the time (or the patience) to manually curate each one to ensure that it was downloaded correctly.

It is good practice to check data integrity when moving files from place to place and there are useful functions for checking data integrity. The two main methods that are used are shasum and md5. shasum is perhaps a bit more common as it is available as standard on MacOS. Both methods generate what is called a “checksum”, a hexadecimal string (e.g. “d241941bac307bd853fd21945d029e62c83cea71”) that is unique to a given file.

The data that you obtained for Topic 1 also contained a file called ```SalmonData_checksums.sha```. If you inspect the contents of ```SalmonData_checksums.sha``` you’ll notice that there are two columns in the file, the left hand column contains a bunch of checksums, the righthand column contains the name of corresponding files. You can compare all the files in the directory you downloaded using:

Store your checksums with your data, whereever you host it.

Navigate to a location on the
```
shasum -c SalmonData_checksums.sha
```
