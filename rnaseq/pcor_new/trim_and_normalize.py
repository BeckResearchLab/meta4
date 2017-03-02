import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

import itertools  # for color palette cycling
import numpy as np
import os
import re
import sys
from cycler import cycler

import pandas as pd
import seaborn as sns

from sklearn.preprocessing import StandardScaler
from sklearn import covariance


RAW_DATA_PATH = '/work/rnaseq/alignments/map_to_contigs_longer_than_1500bp/map_to_contigs_longer_than_1500bp.tsv'
FASTQ_READ_COUNTS_PATH = '/work/rnaseq/rna_fastq_read_counts.tsv' 
PLOT_DIR = './figures'

def load_and_drop_rare_features(min_percent):
    # Load raw counts
    counts = pd.read_csv(RAW_DATA_PATH, sep='\t')
    # Remove underscor columns like __no_feature
    underscores = [c for c in counts['locus'] if '__' in c]
    genes_before = counts.shape[0]
    counts = counts[ ~ counts['locus'].isin(underscores)]
    print('removed {} rows corresponding to htseq-count __ columns like __no_feature'.format(genes_before - counts.shape[0]))

    fastq_counts = pd.read_csv(FASTQ_READ_COUNTS_PATH, sep='\t', names=['sample', 'fastq reads'])
    # remove the.fastq.gz suffix:
    fastq_counts['sample'] = fastq_counts['sample'].str.rstrip('.fastq.gz')

    # Divide by total number of reads in fastq
    cT = counts.set_index('locus').T.reset_index()
    # now columns are genes and rows are samples
    rows_before = cT.shape[0]
    cT.rename(columns={'index':'sample'}, inplace=True)
    cT = pd.merge(cT, fastq_counts, how='outer')
    # divide by the new 'fastq reads' column
    assert cT.shape[0] == rows_before
    # didn't work: cT.div(cT['fastq reads'], axis=0)
    cT.set_index('sample', inplace=True)
    cT = cT.div(cT['fastq reads'], axis=0)
    cT = cT*100  # make it percentages
    # if it worked, there should only be one unique value
    assert cT['fastq reads'].drop_duplicates().shape[0] == 1, 'problem dividing by fastq reads'
    assert cT['fastq reads'].drop_duplicates().iloc[0] == 100, 'problem dividing by fastq reads'
    del cT['fastq reads']

    counts = cT.T  # put the genes back as rows, not columns
    # `counts.sum(axis=1).max()`
    sums = counts.sum(axis=1)
    plot_distribution(sums)

    # trim out genes that don't represent at least x% of the reads in at least one sample 
    shape_before = counts.shape
    genes_to_keep = counts[sums > min_percent]
    print('trimmed genes to set with at least {}% of reads in at least one sample: '
          'shape {} --> {}'.format(min_percent, shape_before, counts.shape))
    genes_to_keep.to_csv('top_features_percent_of_fastq_reads.tsv', sep='\t')
    return genes_to_keep

def plot_distribution(series):
    """
    plot the distriution of sum of percentages
    """
    if not os.path.exists(PLOT_DIR):
        os.mkdir(PLOT_DIR)
    fig, ax = plt.subplots(1, 1, figsize=(4, 2.5))
    series.plot.hist(bins=100, ax=ax)
    figpath = os.path.join(PLOT_DIR, '170302_distiribution_of_percentages_attracted_to_genes.pdf')
    ax.set_yscale('log')
    ax.set_xlabel('sum(percent abundance) across all samples')
    ax.set_ylabel('number of genes')
    fig.savefig(figpath, bbox_inches='tight')

def normalize(counts):
    """
    normalize a counts matrix with samples as rows and genes as columns
    """
    print('normalize the data')
    assert counts.shape[0] < counts.shape[1], 'row samples, and column features expected.  shape: {}'.format(counts.shape)
    print('apply StandardScalar')
    ss = StandardScaler(with_mean=False, with_std=True) 
    scaled = ss.fit_transform(counts)  #shape [n_samples, n_features] 
    # now it is a numpy array. 
    # TODO: write out the column names to re-associate the results with features
    scaled_df = pd.DataFrame(counts, columns=counts.columns, index=counts.index)
    scaled_df.to_csv('top_features_scaled--column_features.tsv', sep='\t')
    print('min value after scaling: {}'.format(scaled_df.min(axis=0).min()))
    print('max value after scaling: {}'.format(scaled_df.max(axis=0).max()))
    return scaled

def graph_lasso(scaled_features, alphas=10.0**np.arange(-3,0)):
    models=dict()
    for alpha in alphas:
        try:
            gl = covariance.GraphLasso(assume_centered=False, verbose=True)
            gl.fit(scaled_features)
        except FloatingPointError:
            print("Failed at alpha = %s" % alpha)
    return models
    



if __name__ == '__main__':
    assert sys.version_info >= (3,0), 'this script expects python 3 to be used'
    
    important_gene_percents = load_and_drop_rare_features(min_percent=0.1)
    genes_scaled_numpy = normalize(important_gene_percents.T)
    #graph_lasso(genes_scaled_numpy)



