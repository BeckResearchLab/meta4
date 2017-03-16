library(DESeq2)

sampleInfoPath = '/work/m4b_binning/assembly/data/sample_info/sample_info_w_cryptic.tsv'
#[ec2-user@ip-10-0-0-158 sample_info]$ head -n 3 sample_info_w_cryptic.tsv
#sample id       LakWas type name        sample number   cryptic metagenome name cryptic metatranscriptome name  oxygen  replicate       week    project
#100_LOW12       LakWasM100_LOW12        100     8777.2.112196.ATCCTA    8903.5.115182.ATCCTA    low     4       12      1056214

# Load in the un-normalized table with all samples and all organisms.  
tsvFile <- '/work/rnaseq/alignments/map_to_contigs_longer_than_1500bp/map_to_contigs_longer_than_1500bp.tsv'
# with locus_tag and product: '/work/rnaseq/alignments/map_to_contigs_longer_than_1500bp/map_to_contigs_longer_than_1500bp_genes.tsv'
# one without a product column:
# /work/rnaseq/alignments/map_to_contigs_longer_than_1500bp/map_to_contigs_longer_than_1500bp.tsv 

# On Waffle, first and 2nd col names were genome	locus_tag	product
# E.g. Acidovora-69x (UID4105) 	Ga0081644_10011	serine/threonine protein kinase
# Row names have to be unique
masterD <- read.table(tsvFile, sep="\t", header=T, quote="", row.names=1)
sampleInfo <- read.table(sampleInfoPath, sep="\t", header=T, quote="", row.names=4)

countData <- masterD
# Remove __ columns
print('Remove __ type columns: ')
print(tail(rownames(masterD), 10))
countData <- countData[ ! grepl('__', rownames(masterD)), ]


head(countData, 2)
# Delete the genome and product columns; they aren't read counts and DESeq doesn't want them. 
#countData$genome <- NULL
#countData$product <- NULL

# remove all 0 rows
# meta4: dim(countData) --> [1] 921440     88
# dim(countData[rowSums(countData) > 0, ]) --> [1] 754840     88
countData <- countData[rowSums(countData) > 0, ]

countData[0:5, 0:5]

# Read in the experiment info; imperative for DESeq's normalization scheme. 
colData <-  sampleInfo[, c('oxygen', 'week')]
head(colData)
# Need to convert 'week' into a factor. 
# d$a <- factor(d$a)
colData$week <- factor(colData$week)

dds <- DESeqDataSetFromMatrix(countData = countData,
                              colData = colData,
                              design = ~ week + oxygen)
vsd <- varianceStabilizingTransformation(dds)
head(assay(vsd), 3)
normCounts <- assay(vsd)

write.table(normCounts, file='normalized_counts.tsv', quote=FALSE, sep='\t')
