library("optparse") # command line arguments
library(ggplot2) # plotting at the end
library(scales) # plotting at the end
library(GeneNet)
rm(list=ls(all=TRUE))

library(reshape)

option_list = list(
  make_option(c("-f", "--file"), type="character", dest="file",
    help="input counts file"),
  # make_option(c("-p", "--percent"), type="float", dest="percent",
  #   help="gene inclusion threshod: minimum percent of reads in at least one sample")
  make_option(c("-p", "--percent"), type="character", dest="percent",
    help="gene inclusion threshod: minimum percent of reads in at least one sample")
);

print('run OptionParser')
opt_parser = OptionParser(option_list=option_list);
print('run parse_args')
opt = parse_args(opt_parser);

print(opt$file)
print(opt$percent)
print(as.numeric(opt$percent))

print('check for file name having been supplied')
if (is.null(opt$file)) {
  print_help(opt_parser)
  stop("Need to supply the input file name, e.g. input_for_R--min_percent_0.005--unnormalized.tsv", call.=FALSE)
}

INPUT = opt$file
PERCENT = as.numeric(opt$percent) # note: the loaded data is already in percents.  No need to divide by 100 again.

print(paste('omit genes with less than', PERCENT, 'perent of reads in at least one sample'))

dir.create('./results', showWarnings = FALSE)
DATA_SAVE_DIR = './results/data/'
dir.create(DATA_SAVE_DIR, showWarnings = FALSE)

GENE_NAMES_TSV='/work/m4b_binning/assembly/prokka/contigs/contigs_longer_than_1500bp/contigs_longer_than_1500bp_gffs_concatenated.gff.genes.tsv'
print(paste('read in gene names:', GENE_NAMES_TSV))
gene_names <- read.csv(GENE_NAMES_TSV, sep='\t', stringsAsFactors=FALSE)

print(paste('read in gene counts:', INPUT))
df <- read.csv(INPUT, sep='\t', stringsAsFactors=FALSE)

df <- rename(df, c(X="locus"))
rownames(df) <- df$locus
# now drop the df column.

loci <- df$locus
df$locus <- NULL

print('transpose dataframe prior to fitting')
df <- t(df)
dim(df)

dfm <- as.matrix(df)
print(dfm[0:4,0:4])
print('preview of colSums:')
print(colSums(dfm)[0:5])
keep_cols <- apply(df, 2, max) > PERCENT
print("keep_cols[0:5]:")
print(keep_cols[0:5])

# trim columns now
print('shape before trimming columns:')
print(dim(dfm))
dfm <- dfm[, keep_cols]
print('shape after trimming columns:')
print(dim(dfm))

# estimate partial correlations
df.static <- ggm.estimate.pcor(dfm, method = "static")  # not time series
# clearly called `estimate.lambda` in corpcor package, but not evident from function signatures

print('test edges')
df.edges <- network.test.edges(df.static, direct=TRUE)
# took 138GB of memory to store df.edges for 26k nodes.

# The memory requirement drops down to 35GB after this extract.network call.  Is it changing the df.edges object?
NUM_EDGES_KEEP = 1e6
print(paste('get top', NUM_EDGES_KEEP, 'edges'))
df.net <- extract.network(df.edges, method.ggm="number", cutoff.ggm=NUM_EDGES_KEEP)
# https://cran.r-project.org/web/packages/GeneNet/GeneNet.pdf
# cutoff.ggm: default cutoff for significant partial correlations

print(paste('shape of df:', dim(df.net)))

node_df <- data.frame(locus=loci, node=c(1:length(loci)), stringsAsFactors=FALSE)
head(node_df)

print('merge on gene names')
# Merge on the locus numbers (e.g. 1_00666) to the network.
node1_names <- setNames(node_df,  c("node1_locus", "node1"))
node2_names <- setNames(node_df,  c("node2_locus", "node2"))
df.net <- merge(df.net, node1_names, all.x=TRUE)
df.net <- merge(df.net, node2_names, all.x=TRUE)

# merge on the gene product names
# prep the shortened locus name for merging
print('merge on locus names')
gene_names$locus <- gsub('contigs_longer_than_1500bp_group_','',gene_names$ID)
gene_names_node1 <- data.frame(node1_locus=gene_names$locus, product_1=gene_names$product)
gene_names_node2 <- data.frame(node2_locus=gene_names$locus, product_2=gene_names$product)

df.net <- merge(df.net, gene_names_node1, all.x=TRUE)
df.net <- merge(df.net, gene_names_node2, all.x=TRUE)
print(head(df.net))

# Save file with top edges and gene names merged on as new tsv
fname <- paste0(DATA_SAVE_DIR, gsub('.tsv', '', INPUT), paste0('--top_', NUM_EDGES_KEEP, '_edges.tsv'))
print(paste('save results for', INPUT, ' (top', NUM_EDGES_KEEP, ' edges) to', fname))
write.table(df.net, file=fname, quote=FALSE, sep='\t')

#-------- Done with real computation.  Now gather some fun facts ----
print('Done with real computation.  Now gather some fun facts')
num_nodes_and_edges_in_trimmed_network <- function(network_df){
	nodes1 <- unique(network_df$node1_locus)
	nodes2 <- unique(network_df$node2_locus)
	num_nodes <- length(unique(c(nodes1, nodes2)))
	num_edges <- dim(network_df)[1]
	return(list('nodes'=num_nodes, 'edges'=num_edges))
}

num_nodes_and_edges <- function(network_df, num_edges){
	print(paste('--- gather stats for network with', num_edges, 'edges ---'))
	trimmed <- extract.network(network_df, method.ggm="number", cutoff.ggm=num_edges)
	nodes1 <- unique(trimmed$node1)
	nodes2 <- unique(trimmed$node2)
	num_nodes <- length(unique(c(nodes1, nodes2)))
	num_edges <- dim(trimmed)[1]

	# min max, median, mean:
	pcor_magnitude_median <- median(abs(trimmed$pcor))
	pcor_magnitude_mean <- mean(abs(trimmed$pcor))
	pcor_min <- min(trimmed$pcor)
	pcor_max <- max(trimmed$pcor)

	return(list('nodes'=num_nodes, 'edges'=num_edges,
		'pcor_magnitude_median'=pcor_magnitude_median, 'pcor_magnitude_mean'=pcor_magnitude_mean,
		'pcor_min'=pcor_min, 'pcor_max'=pcor_max))
	}

# do.call only accepts one arg.  So make something like a Python partial.
num_nodes_and_edges_df_edges <- function(num_edges){
	return(num_nodes_and_edges(df.edges, num_edges))
}

#n Loop over different numbers of edges and compute the number of nodes,
#sizes <- c(1e2, 1e3)  # for dev mode
# max size of network:
num_edges_possible =  dim(df.edges)[1]
nep_rounded = round(log(num_edges_possible,10), 1)
print(paste0('max size of network: 1e', nep_rounded, '10'))
sizes <- c(1e2, 1e3, 1e4, 1e5, 5e5, 1e6, 5e6, 1e7, 1e8, 1e9, 1e10)
print('sizes that are too big (too many edges for total # possible:')
print(sizes[which(sizes>num_edges_possible)])
sizes <- sizes[ - which(sizes>num_edges_possible)]
print('sizes to try:')
print(sizes)
summary <- data.frame(do.call(rbind, Map(num_nodes_and_edges_df_edges, sizes)))

# they elements are in lists.  Unlist them.
summary <- apply(summary, 2, unlist)
summary <- data.frame(summary)
print('summary of nodes, edges at different cutoffs:')
print(head(summary, 2))

fname <- paste0(DATA_SAVE_DIR, gsub('.tsv', '', INPUT), paste0('--summary_for_different_cutoffs.tsv'))
print(paste('save summary of networks at different cutoffs as ', fname))
write.table(summary, file=fname, quote=FALSE, sep='\t')

# plot results
print('plot with ggplot')
p <- ggplot(summary, aes(x=nodes, y=edges, size=pcor_magnitude_median, color=pcor_magnitude_median))
p <- p + geom_point() +
	scale_y_log10(labels = trans_format("log10", math_format(10^.x)),
			breaks = trans_breaks("log10", function(x) 10^x)) +
	scale_x_log10(labels = trans_format("log10", math_format(10^.x)),
                        breaks = trans_breaks("log10", function(x) 10^x))
	#theme(legend.position="bottom", legend.direction='horizontal', legend.box = "vertical")
#ggsave('./results/summary_nodes_edges.pdf', width=5, height=5)
ggsave('./results/summary_nodes_edges.pdf', width=6, height=4)
