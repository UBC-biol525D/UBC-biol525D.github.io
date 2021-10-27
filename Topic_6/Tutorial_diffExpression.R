rm(list=ls())

## Commented out because I've already installed the program

#source("http://bioconductor.org/biocLite.R")
#BiocManager::install("DESeq2")

library("DESeq2")

# Set working directory - This will be different on your machine
directory <- "~/UBC/ZoologyStuff/BIOL525D/Topic_X_RNA/read_counts/"
setwd(directory)

# Set the prefix for each output file name
outputPrefix <- "RNA_tutorial"

# List the names for each of the count files
sampleFiles<- c("warm_sample_01.read_counts.c.txt",
                "warm_sample_03.read_counts.c.txt",
                "warm_sample_02.read_counts.c.txt",
                "warm_sample_04.read_counts.c.txt",
                "warm_sample_05.read_counts.c.txt",
                "warm_sample_06.read_counts.c.txt",
                "cold_sample_07.read_counts.c.txt",
                "cold_sample_08.read_counts.c.txt",
                "cold_sample_09.read_counts.c.txt",
                "cold_sample_10.read_counts.c.txt",
                "cold_sample_11.read_counts.c.txt",
                "cold_sample_12.read_counts.c.txt")

# List the sample names
sampleNames<- c("sample_01",
                "sample_02",
                "sample_03",
                "sample_04",
                "sample_05",
                "sample_06",
                "sample_07",
                "sample_08",
                "sample_09",
                "sample_10",
                "sample_11",
                "sample_12")

# A list of treatments
sampleCondition <- c(rep("warm", 6),rep("cold", 6))

# Combine the sample info into a table for DESeq
sampleTable <- data.frame(sampleName = sampleNames, fileName = sampleFiles, condition = sampleCondition)

# A list of the factors for a future dataframe
treatments = c("warm","cold")



# Built the DESeq dataset object using the counts data
ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                       directory = directory,
                                       design = ~ condition)

# Estimate size factors - normalize read counts
dds <- estimateSizeFactors(ddsHTSeq)

# Take a look at the raw counts
View(counts(ddsHTSeq))
# Take a look at the normalised counts
View(counts(dds, normalized=TRUE))

##########################################

# Convert treatment into factors
colData(ddsHTSeq)$condition <- factor(colData(ddsHTSeq)$condition,
                                      levels = treatments)


# The guts of the analysis - run the DESeq2 GLM on the count data
dds <- DESeq(ddsHTSeq)

res <- results(dds)
# Make a simple volcano plot
plot(res$log2FoldChange, -log10(res$padj))
library(ggplot2)

res_DF <- data.frame(res)

# A quick Volcano plot
ggplot(data = res_DF, aes( y = -log10(padj), x = log2FoldChange, fill = abs(log2FoldChange)))+
  geom_point(shape = 16)+
  geom_point(shape = 21)+
  scale_fill_gradient(low="white",high = "red")+
  theme_bw()

# A quick Volcano plot - log10 the y-axis
ggplot(data = res_DF, aes( y = -log10(padj), x = log2FoldChange, fill = abs(log2FoldChange)))+
  geom_point(shape = 16)+
  geom_point(shape = 21)+
  scale_fill_gradient(low="white",high = "red")+
  scale_y_log10()+
  theme_bw()



# transform raw counts into normalized values
# DESeq2 has two options:  1) rlog transformed and 2) variance stabilization
# variance stabilization is very good for heatmaps, etc.
rld <- rlogTransformation(dds, blind=T)
vsd <- varianceStabilizingTransformation(dds, blind=T)

#Principal components plot shows additional but rough clustering of samples
library("genefilter")
library("ggplot2")
library("grDevices")

rv <- rowVars(assay(rld))
select <- order(rv, decreasing=T)[seq_len(min(500,length(rv)))]
pc <- prcomp(t(assay(vsd)[select,]))


# set condition
condition <- treatments
scores <- data.frame(pc$x)
scores$condition <- rep(c("warm","cold"),each =6)
(pcaplot <- ggplot(scores, aes(x = PC1, y = PC2, col = (factor(condition))))
  + geom_point(size = 5)
  + ggtitle("Principal Components")
  + scale_colour_brewer(name = " ", palette = "Set1")
  + theme(
    plot.title = element_text(face = 'bold'),
    legend.key = element_rect(fill = 'NA'),
    legend.text = element_text(size = 10, face = "bold"),
    axis.text.y = element_text(colour = "Black"),
    axis.text.x = element_text(colour = "Black"),
    axis.title.x = element_text(face = "bold"),
    axis.title.y = element_text(face = 'bold'),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.background = element_rect(color = 'black',fill = NA)
  ))

ggsave(pcaplot,file=paste0(outputPrefix, "-ggplot2.pdf"))

