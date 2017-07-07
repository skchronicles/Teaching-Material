# Example script for Bioinfo class
# First -- change working directory to desired folder containing .Cel and phenotype files!
# if needed, download required tools by these commands:
source("http://www.bioconductor.org/biocLite.R")
biocLite("affy") 
biocLite("affycoretools") 
library(affycoretools)

#Loading affy and affycoretools packages into our R environment: 
library(affy) 
library(affycoretools)

# read .cel files into "mydata" affy object
mydata <- ReadAffy()

#read phenotype (sample identify) file
pData(mydata)<-read.table("PhenoData_FTY.txt", header=T, row.names=1, sep="\t")

# check content of ”pdata”
pData(mydata)

#check content of “mydata”
mydata

#check distribution of data, quality of files
hist(mydata)
boxplot(mydata)
image(mydata [,1])

# run quantile normalization, RMA summarization, and name file "eset"
eset <- rma(mydata)
# check histogram/boxplot of RMA output
hist(exprs(eset))
boxplot(exprs(eset))

# write RMA results to a tab delimited file called ”ftytest.txt"
write.exprs(eset,file="ftytest.txt", sep="\t")

#Alternative approach for outputting data
#RMAresults=exprs(eset)
#write.table(RMAresults,file="ftytest.txt", sep="\t")

#plot PCA of data to check for batch effects
#This method is not working at present -- need to debug, likely deprecated
#biocLite("BioCGenerics")
#plotPCA(eset, groups = as.numeric(pData(mydata1) [,1]), groupnames = levels(pData(mydata1) [,1]))