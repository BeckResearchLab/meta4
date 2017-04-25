# Import matplotlib before seaborn
import matplotlib as mpl
import matplotlib.pyplot as plt

import re
import sys

import numpy as np
import pandas as pd

from scipy.stats import mannwhitneyu

# for a regex:
sys.path.append('/work/general_scripts')
from plot_subplots import CONTIG_PARSE_REGEX

def load_edges(edge_file):
    edges = pd.read_csv(edge_file, sep='\t')
    return edges

def extract_second_num(string):
    """
    5_151185 --> 151185  (int type)
    To know whether two genes are adjacent or close on a contig, we need
    to compare the gene number.
    E.g. 151185 and 151186 should be adjacent.
    """
    m = re.search('[0-9]+_([0-9]+)', string)
    if m:
        found = m.group(1)
        return int(found)

def round_to_nearest_step(x, step_width):
    """
    For histograms.  E.g. 0.032 --> 0.03 for inputs 0.032, 0.01
    """
    return round(float(x) / step_width) * step_width

def make_bin_edges(min_val, max_val, bin_width=0.02):
    """
    The histogram edges have one boundary at zero so + and - signed
    pcors are kept separate.
    """
    lower_boundary = round_to_nearest_step(min_val, bin_width)
    upper_boundary = round_to_nearest_step(max_val, bin_width)
    print('lower: {}, upper: {}'.format(lower_boundary, upper_boundary))
    return np.arange(lower_boundary, upper_boundary, bin_width)

def get_edges(edge_df, edge_string):
    selected_edge_df = \
    edge_df[edge_df['product_1'].str.contains(edge_string) &
          edge_df['product_2'].str.contains(edge_string)]
    num_selected_edge_df = selected_edge_df.shape[0]
    print('number of edges selected: {}'.format(num_selected_edge_df))
    return selected_edge_df

def get_mmo_edges(edge_df):
    return get_edges(edge_df, 'ethane mono')

def get_h6p_edges(edge_df):
    return get_edges(edge_df, '3-hexulose-6-phosphate')

def get_mdh_subunit_edges(edge_df):
    return get_edges(edge_df, 'Methanol dehydrogenase')

def separate_out_linked_and_not_linked_rows(edge_df, distance_allowed):
    """
    Approximate which genes "belong" together (e.g. are in an operon)
    based on being with the specified distance.
    Store the pcor values in each category in a list.
    """
    linked = []
    unlinked = []

    for index, row in edge_df.sort_values(['node1_locus', 'node2_locus']).iterrows():
        # TODO: make sure they are also in the same batch. (E.g. 1_100 and
        # 2_101 should not be considered adjacent.)
        # TODO: make sure they are on the same contig.  It's possible
        # 1_1000 and 1_1001 can be on different contigs.
        n1 = extract_second_num(row['node1_locus'])
        n2 = extract_second_num(row['node2_locus'])
        if abs(n1-n2) <= distance_allowed:
            print(row['node1_locus'], row['node2_locus'], row['pcor'])
            print('    {} | {}'.format(row['product_1'], row['product_2']))
            linked.append(row['pcor'])
        else:
            unlinked.append(row['pcor'])
    print('Mann–Whitney U test: {}'.format(test_mannwhitneyu(linked, unlinked)))
    return linked, unlinked

def test_mannwhitneyu(x, y):
    """
    Run the scipy stats implementation.
    https://en.wikipedia.org/wiki/Mann%E2%80%93Whitney_U_test
    """
    return mannwhitneyu(x, y)

def filter_out_same_same_edges(edge_df):
    """
    Remove rows in the data frame if the two gene products are identical.
    E.g. remove out hps:hps pairs if you are looking for hps:hpi pairs.
    """
    for index, row in edge_df.sort_values(['node1_locus', 'node2_locus']).iterrows():
        p1 = row['product_1']
        p2 = row['product_2']
        if p1 == p2:
            print('drop row with products "{}":"{}"'.format(p1, p2))
            edge_df.drop(index, inplace=True)

def filter_out_edges_with_same_node(edge_df, node_strings):
    # first get the
    for ns in node_strings:
        for index, row in edge_df.sort_values(['node1_locus', 'node2_locus']).iterrows():
            p1 = row['product_1']
            p2 = row['product_2']
            if (ns in p1) and (ns in p2):
                print('drop row with products "{}":"{}"'.format(p1, p2))
                edge_df.drop(index, inplace=True)
    # no need to return b/c deleted in place.

def plot_linked_and_unlinked(linked, unlinked, orientation='vertical',
                             bin_label='partial correlation', bin_width=0.015 ):
    """
    Make a stacked histogram from the two lists of pcor values.
    """
    min_val = min(min(linked), min(unlinked))
    max_val = max(max(linked), max(unlinked))
    print(min_val, max_val)
    bins = make_bin_edges(min_val, max_val, bin_width)
    print('bin edges: {}'.format(bins))

    if orientation == 'vertical':
        figsize = (3.5, 5)
    else:
        figsize = (6, 3)

    fig, ax = plt.subplots(1,1, figsize=figsize)
    ax.hist([linked, unlinked],
             bins=bins, stacked=True, color = ['#31a354', '#bdbdbd'],
             label=['operon pair', 'non-operon pair'], orientation=orientation)
    plt.xticks(rotation=90)
    if orientation == 'vertical':
        plt.axvline(0, color='black', linestyle='-', linewidth=.5)
        ax.set_xlabel(bin_label)
        ax.set_ylabel('count')
    else:
        plt.axhline(0, color='black', linestyle='-', linewidth=.5)
        ax.set_ylabel(bin_label)
        ax.set_xlabel('count')

    plt.legend()
    return fig


#-----------  Now for asking more generally whether same-contig gene pairs have higher pcor values.

def get_genes_contigs():
    GENES_CONTIG_FILE = '/work/m4b_binning/assembly/prokka/contigs/contigs_longer_than_1500bp/contigs_longer_than_1500bp_gffs_concatenated.gff.genes.tsv'
    genes = pd.read_csv(GENES_CONTIG_FILE, sep='\t')
    genes['locus'] = genes['ID'].str.extract(CONTIG_PARSE_REGEX, expand=True)
    return genes

def merge_contig_info_onto_edge_df(edges):
    genes = get_genes_contigs()
    cd1 = {'contig':'node1_contig', 'locus':'node1_locus'}
    cd2 = {'contig':'node2_contig', 'locus':'node2_locus'}
    edges = pd.merge(edges,
                      genes[['contig', 'locus']].rename(columns=cd1))
    edges = pd.merge(edges,
                      genes[['contig', 'locus']].rename(columns=cd2))
    edges['same contig'] = edges['node1_contig'] == edges['node2_contig']
    return edges

def plot_pdf_of_same_and_different_contig_pcors(edge_df):
    fig, ax = plt.subplots(1, 1, figsize=(4, 2.5),
                        sharey=False, sharex=True)
    same_contig_pcors = edge_df[edge_df['same contig']]

    same = edge_df[edge_df['same contig'] == True]
    print("same.shape: {}".format(same.shape))
    not_same = edge_df[edge_df['same contig'] == False]
    assert edge_df.shape[0] == same_contig_pcors.shape[0] + not_same.shape[0]
    print("not_same.shape: {}".format(not_same.shape))
    test_mannwhitneyu(same['pcor'], not_same['pcor'])

    a=0.3
    #bins = np.arange(-0.04, 0.14, 0.004)
    bins = np.arange(-0.01, 0.03, 0.001)
    #bins = np.arange(data1.min(), data2.max(), 2000)
    ax.hist(same['pcor'], bins, alpha=a, color='green', normed=True)
    ax.hist(not_same['pcor'], bins, alpha=a, color='blue', normed=True)
    ax.set_xlabel('partial correlation (GeneNet)')
    ax.set_ylabel('pdf (normed count)')
    plt.xticks(rotation=90)
    #ax.set_yscale('log') # looks funny.

    return fig

