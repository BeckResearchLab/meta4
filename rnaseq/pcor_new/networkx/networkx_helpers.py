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

def add_edges_from_df(network, df):
    startTime = datetime.now()
    for idx, row in df.iterrows():
    # add_edge: The nodes u and v will be automatically added if
    # they are not already in the graph.
        network.add_edge(row['gene A'], row['gene B'],
                    attr_dict = {'pcor': row['pcor']} )
    print('networkx edge adding time: {}'.format(datetime.now() - startTime))
    return network

def get_unique_nodes(df):
    df_nodes = df[['gene A']].drop_duplicates()
    df_nodes.rename(columns={'gene A':'gene'}, inplace=True)
    df_nodes = \
        df_nodes.merge(df[['gene B']].drop_duplicates().rename(
            columns={'gene B':'gene'}))
    return df_nodes

def build_network(edges_path, tail_percent, genes_path):
    edges = load_edges(edges_path, tail_percent)
    genes = load_gff_tsv(genes_path)
    genes_in_edges = get_unique_nodes(edges)
    genes = genes[genes['ID'].isin(genes_in_edges['gene'])]

    print('number of unique nodes: {}'.format(genes.shape[0]))

    #n = networkx.DiGraph()  # directed
    n = networkx.Graph()   # undirected

    n = add_nodes_from_df(network=n, df=genes)
    n = add_edges_from_df(network=n, df=edges)
    print('number of nodes: {}'.format(len(n.nodes())))
    print('number of edges: {}'.format(len(n.edges())))

    return n
