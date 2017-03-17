import datetime
import numpy as np
import pandas as pd

# Import matplotlib before seaborn
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns

# dev-scale:
#fname = 'Ledoit_Wolf_r4.8xlarge_244GB_memory/ledoit_wolf_precision_cutoff_0.05--4465_genes.tsv'
# full scale
fname = 'Ledoit_Wolf_r4.8xlarge_244GB_memory/ledoit_wolf_precision_cutoff_0.005--26967_genes.tsv'

gene_names_path = '/work/m4b_binning/assembly/prokka/contigs/contigs_longer_than_1500bp/contigs_longer_than_1500bp_gffs_concatenated.gff.genes.tsv'

all_edges = pd.read_csv(fname, sep='\t')

print('sort by pcor')
all_edges.sort_values('pcor', ascending=False, inplace=True)
print('done sorting')

#percentage = 2.5/100.
percentage = 5
quantile_tails = percentage/100/2
extremes = all_edges[(all_edges['pcor'] >= all_edges['pcor'].quantile(1 - quantile_tails)) |
              (all_edges['pcor'] <= all_edges['pcor'].quantile(quantile_tails))]

rows, cols = extremes.shape
#nodes = extremes[['gene A','gene B']].apply(np.unique)
print('calculate # of nodes')
nodes = len(np.unique(all_edges[['gene A','gene B']]))
print('# of nodes is calculated')

# prep output file name
result_fname = fname.replace('.tsv', '--extreme_{}_percent--{}_nodes_{}_edges.tsv'.format(percentage, nodes, rows))
print('save results to', result_fname)

# merge on the gene info.
gene_names = pd.read_csv(gene_names_path, sep='\t')
# contigs_longer_than_1500bp_group_1_00001 --> 1_00001
gene_names['ID'] = gene_names['ID'].str.extract('[A-z0-9_]+([0-9]+_[0-9]+)', expand=True)
gene_names_A = gene_names.copy().rename(columns={'ID': 'gene A', 'product':'product A', 'bp':'bp A'})
gene_names_B = gene_names.copy().rename(columns={'ID': 'gene B', 'product':'product B', 'bp':'bp B'})
extremes = pd.merge(extremes, gene_names_A, how='left')
extremes = pd.merge(extremes, gene_names_B, how='left')
assert extremes.shape[0] == rows, 'uh-oh: merge problem'

extremes.to_csv(result_fname, sep='\t', index=False)

# Make a hist of the original and new.
print('plot the histograms of all edges, and the selected extreme edges')
fig, ax = plt.subplots(1, 1, figsize=(4,3))
all_edges['pcor'].rename('un-trimmed').hist(ax=ax, bins=100, color='gray')
extremes['pcor'].rename('trimmed').hist(ax=ax, bins=100, color='#74c476')
#ax.set_xscale('log') # can't take log of a negative number.
ax.set_yscale('log')
ax.set_xlabel('precision matrix value')
ax.set_ylabel('number of edges')
fig.savefig(result_fname.replace('.tsv', '.pdf'), bbox_inches='tight')
print('done plotting')
