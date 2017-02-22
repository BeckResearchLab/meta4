library(DESeq2)

# Load in the un-normalized table with all samples and all organisms.  
tsvFile <- '../assemble_summaries/summary.dat'

print(tsvFile)

masterD <- read.table(tsvFile, sep="\t", header=T, quote="", row.names=2)

countData <- masterD
head(countData)
# Delete the genome and product columns; they aren't read counts and DESeq doesn't want them. 
countData$genome <- NULL
countData$product <- NULL

# remove all 0 rows
# Was: (160226)
# countData <- countData[rowSums(countData[, -1]) > 0, ]

# the colnames at this point: 
# > colnames(countData)
#  [1] "LakWasM100_LOW12_2_rpkm" "LakWasMe97_LOW12_2_rpkm"
#  [3] "LakWasMe98_LOW12_2_rpkm" "LakWasMe99_LOW12_2_rpkm"
#  [5] "LakWasM104_HOW12_2_rpkm" "LakWasM105_HOW12_2_rpkm"
# What happens when we do the `countData[, -1]` ?? 
# Removes the first column.  This is NOT The desired behavior if both genome and product are set to NULL.
head(countData)

# Read in the experiment info; imperative for DESeq's normalization scheme. 
# ---------FIX: added sample_info to current dir for development on local computer-----------------
colData <- read.table("../sample_info.xls", sep="\t", header=T, quote="", row.names=1)
#colData <- read.table("./sample_info.xls", sep="\t", header=T, quote="", row.names=1)
head(colData)

dds <- DESeqDataSetFromMatrix(countData = countData,
                              colData = colData,
                              design = ~ week + O2)
# remove rows with all zeros before applying DESeq. 
dds <- dds[ rowSums(counts(dds)) > 1, ]
#dds <- DESeq(dds)
#res <- results(dds)
#rld <- rlog(dds)
#head(assay(rld), 3)
vsd <- varianceStabilizingTransformation(dds)
head(assay(vsd), 3)
normCounts <- assay(vsd)

write.table(normCounts, file='normalized_counts.tsv', quote=FALSE, sep='\t')
