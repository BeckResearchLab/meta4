import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

from datetime import datetime
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

import argparse


RAW_DATA_PATH = '/work/rnaseq/alignments/map_to_contigs_longer_than_1500bp/map_to_contigs_longer_than_1500bp.tsv'
FASTQ_READ_COUNTS_PATH = '/work/rnaseq/rna_fastq_read_counts.tsv' 
GENE_NAMES_PATH = '/work/m4b_binning/assembly/prokka/contigs/contigs_longer_than_1500bp/contigs_longer_than_1500bp_gffs_concatenated.gff.genes.tsv'

PLOT_DIR = './figures'

def load_and_drop_rare_features(min_percent, plot_dist=True):
    """
    Load original counts tsv, divide counts by number of reads in RNA-seq fastq,
    delete genes that aren't above min_percent in at least one sample,
    and return it.
    """
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
    # after dividing by the number of fastq reads, and multiplying by 100, the 'fastq reads' 
    # column should now be 100 at every value. 
    assert cT['fastq reads'].drop_duplicates().shape[0] == 1, 'problem dividing by fastq reads'
    assert cT['fastq reads'].drop_duplicates().iloc[0] == 100, 'problem dividing by fastq reads'
    # delete the column of 100s now. 
    del cT['fastq reads']

    counts = cT.T  # put the genes back as rows, not columns
    # `counts.sum(axis=1).max()`
    sums = counts.sum(axis=1)
    if plot_dist:
        plot_distribution(sums)

    # trim out genes that don't represent at least x% of the reads in at least one sample 
    shape_before = counts.shape
    genes_to_keep = counts[sums > min_percent]
    print('trimmed genes to set with at least {}% of reads in at least one sample: '
          'shape {} --> {}'.format(min_percent, shape_before, genes_to_keep.shape))
    return genes_to_keep  # dataframe

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

def shuffle_cols(df):
    """
    Shuffle the columns of a dataframe, before cross-val sub-sampling
    """
    return df.sample(frac=1)
    

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
    # Put the row and column names back on.
    scaled_df = pd.DataFrame(counts, columns=counts.columns, index=counts.index)
    print('min value after scaling: {}'.format(scaled_df.min(axis=0).min()))
    print('max value after scaling: {}'.format(scaled_df.max(axis=0).max()))
    return scaled_df

def ledoit_wolf(scaled_features):
    # sklearn.covariance.LedoitWolf
    # class sklearn.covariance.LedoitWolf(store_precision=True, assume_centered=False, block_size=1000)[source]
    # Ledoit-Wolf optimal shrinkage coefficient estimate
    print("computing Ledoit-Wolf optimal shrinkage coeffecient estimate")
    start_time = datetime.now()
    print('time: {}'.format(str(start_time)))
    lwe = covariance.LedoitWolf().fit(scaled_features.as_matrix())
    end_time = datetime.now()
    total_time = end_time - start_time
    print('LedoitWolf time for {} genes: {}'.format(scaled_features.shape[1], str(total_time)))
    return lwe

def run_ledoit_wolf(genes_scaled, include_product_name_cols=False):
    print(genes_scaled.shape)
    lw = ledoit_wolf(genes_scaled)
    # precision matrix is proportional to the partial correlation matrix. 
    prec = lw.get_precision()
    # These are precision values, not partial correlations (not in range [-1, 1]) 
    print(prec.shape)
    prec_df = pd.DataFrame(prec, 
                           columns=genes_scaled.columns, 
                           index=genes_scaled.columns)
    # http://stackoverflow.com/questions/34417685/melt-the-upper-triangular-matrix-of-a-pandas-dataframe
    # upper_triangle_indices: np.triu(np.ones(prec.shape), 1).astype(np.bool)
    # Set the lower triangle and diagonal to NaN
    prec_triangle = prec_df.where(np.triu(np.ones(prec_df.shape), 1).astype(np.bool))
    prec_triangle = prec_triangle.stack().reset_index()
    prec_triangle.columns = ['gene A', 'gene B', 'pcor']

    # merge on gene products
    print(prec_triangle.head(2))
    if include_product_name_cols:
        print('merge on the gene names to the flattened triangle of precision matrix values.  Storage memory intensive!!')
        print('show memory before merging on colnames:')
        print(prec_triangle.info())
        gene_names = pd.read_csv(GENE_NAMES_PATH, sep='\t')
        for c in ['A', 'B']:
            gn = gene_names.copy()
            gn.rename(columns={'ID':'gene {}'.format(c), 
                               'product': 'product {}'.format(c)}, inplace=True)
            prec_triangle = pd.merge(prec_triangle, gn)
        print(prec_triangle.info())
    else: 
        print('not merging on gene product names, for the sake of memory.  See {}'.format(GENE_NAMES_PATH))
        print('memory used:')
        print(prec_triangle.info())
    return prec_triangle

def ledoit_wolf_with_increasing_size(abundance_cutoff_list, dirname='ledoit_wolf_precision', include_product_name_cols=False):
    if not os.path.exists(dirname):
        os.mkdir(dirname)
    # Loop through the selected values. 
    for ac in abundance_cutoff_list:
        print('start ledoit_wolf for genes with % abundance in at least one sample > {}'.format(ac))
        important_genes = load_and_drop_rare_features(min_percent=ac, plot_dist=False)
        num_genes = important_genes.shape[0]
        normalized = normalize(important_genes.T)
        print('run ledoit wolf for abundance cutoff = {} ({} genes)'.format(ac, num_genes))
        precision = run_ledoit_wolf(normalized, include_product_name_cols=include_product_name_cols) # gets to be a lot of memory!
        filename = 'ledoit_wolf_precision' + '_cutoff_{}--{}_genes.tsv'.format(ac, num_genes)
        filename = os.path.join(dirname, filename)
        print('save results as filename: ', filename)
        precision.to_csv(filename, sep='\t')
        print('----------------------------------------------')
        sys.stdout.flush()


def graph_lasso(scaled_features, alphas=10.0**np.arange(-3,0)):
    models=dict()
    for alpha in alphas:
        print('try graph lasso with alpha = {}'.format(alpha))
        try:
            start_time = datetime.now()
            gl = covariance.GraphLasso(assume_centered=False, verbose=True)
            gl.fit(scaled_features)
            total_time = end_time - start_time
            print('GraphLasso time for matrix with dim {}: {}'.format(
                scaled_features.shape, str(total_time)))
            print('finished alpha = {} without error'.format(alpha))
        except FloatingPointError:
            print("Failed at alpha = %s" % alpha)
    return models
    

if __name__ == '__main__':

    assert sys.version_info >= (3,0), 'this script expects python 3 to be used'

    # get args
    parser = argparse.ArgumentParser()
    parser.add_argument('-l','--list', nargs='+', type=float, help='<Required> list of cutoffs to use', required=True)
    parser.add_argument('-d','--dir_name', help='<Required> directory to store resulting tsv files in', required=True)
    args = parser.parse_args()
    
    #
    print('use cutoff values {}'.format(args.list))
    print('save resultig tsv files to "{}"'.format(args.dir_name))
    assert type(args.list[0]) == float, 'expected float but got {}'.format(type(args.list[0]))

    # 1e-10 was shown to have all 754836 previously.
    #cutoffs = [1, 0.1, 0.05, 0.01, 0.005, 0.001, 0.0005, 0.0001, 5e-4, 1e-4, 5e-5, 1e-5, 1e-10]
    cutoffs = args.list
    print('Run Ledoit Wolf on increasing data sizes, according to cutoffs {}'.format(cutoffs))
    ledoit_wolf_with_increasing_size(cutoffs, dirname=args.dir_name, include_product_name_cols=False)
    
