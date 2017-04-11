# Import matplotlib before seaborn
import matplotlib as mpl
import matplotlib.pyplot as plt
import itertools  # for color palette cycling
import re
import pandas as pd
import seaborn as sns
import sys
from cycler import cycler
import seaborn as sns

from datetime import datetime
import pickle
import timeit

import networkx

GENE_PATH = '/work/m4b_binning/assembly/prokka/contigs/contigs_longer_than_1500bp/contigs_longer_than_1500bp_gffs_concatenated.gff.genes.tsv'


def load_edges(network_path, tail_percent):
    startTime = datetime.now()
    tail_decimal = tail_percent/100.
    df = pd.read_csv(network_path, sep='\t')
    print('time to load un-trimmed edge file: {}'.format(datetime.now() - startTime))
    print('select out the most positive and most negative {} percent of edges'.format(tail_percent))
    extremes = df[(df['pcor'] >= df['pcor'].quantile(1 - tail_decimal)) |
              (df['pcor'] <= df['pcor'].quantile(tail_decimal))]
    print('reduced edges from {} to {}'.format(df.shape[0], extremes.shape[0]))
    return extremes

def load_gff_tsv(tsv_path):
    df = pd.read_csv(tsv_path, sep='\t')
    df['ID'] = df['ID'].str.extract('[A-z0-9]+_([0-9]+_[0-9]+)', expand=True)
    return df

def add_nodes_from_df(network, df):
    # get set of unique nodes.
    startTime = datetime.now()
    # load nodes
    print(df.head())
    for idx, row in df.iterrows():
        network.add_node(n=row['ID'], attr_dict={'product':row['product'], 'contig':row['contig']})
    print('networkx node adding time: {}'.format(datetime.now() - startTime))
    return network

def add_edges_from_df(network, df, gene1colname, gene2colname):
    startTime = datetime.now()
    for idx, row in df.iterrows():
    # add_edge: The nodes u and v will be automatically added if
    # they are not already in the graph.
        network.add_edge(row[gene1colname], row[gene2colname],
                    attr_dict = {'pcor': row['pcor']} )
    print('networkx edge adding time: {}'.format(datetime.now() - startTime))
    return network

def get_unique_nodes(df, gene1colname, gene2colname):
    df_nodes = df[[gene1colname]].drop_duplicates()
    df_nodes.rename(columns={gene1colname:'gene'}, inplace=True)
    df_nodes = \
        df_nodes.merge(df[[gene2colname]].drop_duplicates().rename(
            columns={gene2colname:'gene'}))
    return df_nodes.drop_duplicates()

def plot_edge_dist(pcors):
    fig, ax = plt.subplots(1,1, figsize=(4,3))
    pcors.plot.hist(ax=ax)

def get_loci_colnames(df):
    if 'node1_locus' in df.columns:
        return 'node1_locus', 'node2_locus'
    elif 'gene A' in df.columns:
        return 'gene A', 'gene B'


def build_network(edges_path, tail_percent=None, genes_path=GENE_PATH):
    # TODO: somehow there are a few nodes without gene product attributes.
    # TODO: how did they get in?
    if tail_percent is not None:
        edges = load_edges(edges_path, tail_percent)
    else:
        edges = pd.read_csv(edges_path, sep='\t')

    print(edges.columns)
    gene1colname, gene2colname = get_loci_colnames(edges)

    genes = load_gff_tsv(genes_path)
    genes_in_edges = get_unique_nodes(edges, gene1colname, gene2colname)
    genes = genes[genes['ID'].isin(genes_in_edges['gene'])]

    print('number of unique nodes: {}'.format(genes.shape[0]))

    #n = networkx.DiGraph()  # directed
    n = networkx.Graph()   # undirected

    n = add_nodes_from_df(network=n, df=genes)
    n = add_edges_from_df(n, edges, gene1colname, gene2colname)
    print('number of nodes: {}'.format(len(n.nodes())))
    print('number of edges: {}'.format(len(n.edges())))

    return n

