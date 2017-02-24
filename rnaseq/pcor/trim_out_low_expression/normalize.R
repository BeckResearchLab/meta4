#!/usr/bin/Rscript

library("optparse") # command line args
print('set up options')

option_list = list(
  make_option(c("-s", "--sample_info"), type="character", default=NULL, 
              #help="sample info file with descriptors (e.g. week, O2) corresponding to each volumn of gene_tsv", 
	      metavar="character"),
  make_option(c("--gene_tsv"), type="character", default=NULL, 
              help="tsv file with *one* column uniquely identifying each gene, and one column per sample", 
	      metavar="character")
); 
 
print('Set up OptionParser')
opt_parser = OptionParser(option_list=option_list);
print('parse arguments')
opts = parse_args(opt_parser);

if (is.null(opts$sample_info)){
  print_help(opt_parser)
  print('suggested sample_info: /work/m4b_binning/assembly/data/sample_info/sample_info_w_cryptic.tsv,')
  print('suggested gene_tsv: /work/rnaseq/alignments/map_to_contigs_longer_than_1500bp/map_to_contigs_longer_than_1500bp.tsv')
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
}

# This import dumps a lot of junk to the terminal.
# If the user didn't enter args, they shouldn't see all that.  So be weird and put it later.
library(DESeq2)

sampleInfoPath = opts$sample_info  #'/work/m4b_binning/assembly/data/sample_info/sample_info_w_cryptic.tsv' 
#[ec2-user@ip-10-0-0-158 sample_info]$ head -n 3 sample_info_w_cryptic.tsv
#sample id       LakWas type name        sample number   cryptic metagenome name cryptic metatranscriptome name  oxygen  replicate       week    project
#100_LOW12       LakWasM100_LOW12        100     8777.2.112196.ATCCTA    8903.5.115182.ATCCTA    low     4       12      1056214

# Load in the un-normalized table with all samples and all organisms.  
tsvFile <- opts$gene_tsv  #'/work/rnaseq/alignments/map_to_contigs_longer_than_1500bp/map_to_contigs_longer_than_1500bp.tsv'
# with locus_tag and product: '/work/rnaseq/alignments/map_to_contigs_longer_than_1500bp/map_to_contigs_longer_than_1500bp_genes.tsv'
# one without a product column:
# /work/rnaseq/alignments/map_to_contigs_longer_than_1500bp/map_to_contigs_longer_than_1500bp.tsv 

# On Waffle, first and 2nd col names were genome	locus_tag	product
# E.g. Acidovora-69x (UID4105) 	Ga0081644_10011	serine/threonine protein kinase
# Row names have to be unique
print('read expression tsv')
masterD <- read.table(tsvFile, sep="\t", header=T, quote="", row.names=1)
print('read sample info tsv')
sampleInfo <- read.table(sampleInfoPath, sep="\t", header=T, quote="", row.names=4)

countData <- masterD
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
