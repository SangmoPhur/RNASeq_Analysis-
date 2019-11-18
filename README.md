# RNASeq_Analysis
#rna seq tutorila via website
library(edgeR)
library(limma)
library(Glimma)
library(gplots)
source("https://bioconductor.org/biocLite.R")
source("http://bioconductor.org/biocLite.R")
biocLite("org.Mm.eg.db")
library(RColorBrewer)
#read the data into the R by downloading from geo
seqdata <- read.delim("Downloads/data/3219685/GSE60450_Lactation-GenewiseCounts.txt",stringsAsFactors = FALSE)
#read the sample info into R
sampleinfo <- read.delim("Downloads/data/3219685/SampleInfo.txt")
head(seqdata)
dim(seqdata)
sampleinfo
#remove first two columns from seqdata
countdata <- seqdata[,-(1:2)]
#look at the output
head(countdata)
#store the entrezgeneid as rownames
rownames(countdata) <- seqdata[,1]
head(countdata)
colnames(countdata)
# using substr, you extract the characters starting at position 1 and stopping at position 7 of the colnames
colnames(countdata) <- substr(colnames(countdata),start=1,stop=7)
head(countdata)
table(colnames(countdata)==sampleinfo$SampleName)
#filtering the lowly expressed genes
#cpm counts per million. we keep genes that are expressed at counts per million above 0.5 in at least two samples 
#using the cpm function from edgeR library to generate cpm values and then filter. we are normalizing the different sequencig depth of the samples
#obtaim cpms
myCpm <- cpm(countdata)
#have a look at the output
head(myCpm)
#which values in mycpm are greater than 0.5? 
thresh <- myCpm > 0.5
head(thresh)
#summary of how many trues are there in each row
table(rowSums(thresh))
#we would like to keep genes that have at least 2 trues in each row of the thressh 
keep <- rowSums(thresh) >= 2
#subset the rwos of countdata to keep the more highly expressed genes
keep <- rowSums(thresh) >= 2
#subset the rows of countdata to keep the more hihgly expressed genes
counts.keep <- countdata[keep,]
summary(keep)
dim(counts.keep)
plot(myCpm[,1],countdata[,1])
# Let us limit the x and y-axis so we can actually look to see what is happening at the smaller counts
plot(myCPM[,1],countdata[,1],ylim=c(0,50),xlim=c(0,3))
# Add a vertical line at 0.5 CPM
abline(v=0.5)
######
y <- DGEList(counts.keep)
# have a look at y
y
# See what slots are stored in y
names(y)
# Library size information is stored in the samples slot
y$samples
#quality control
y$samples$lib.size
# The names argument tells the barplot to use the sample names on the x-axis
# The las argument rotates the axis names
barplot(y$samples$lib.size,names=colnames(y),las=2)
# Add a title to the plot
title("Barplot of library sizes")
# Get log2 counts per million
logcounts <- cpm(y,log=TRUE)
# Check distributions of samples using boxplots
boxplot(logcounts, xlab="", ylab="Log2 counts per million",las=2)
# Let's add a blue horizontal line that corresponds to the median logCPM
abline(h=median(logcounts),col="blue")
title("Boxplots of logCPMs (unnormalised)")
