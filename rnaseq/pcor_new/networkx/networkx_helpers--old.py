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

def load_df(fname):
    startTime = datetime.now()
    dataframe = pd.read_csv(fname, sep='\t')
    print('shape of input data: {}'.format(dataframe.shape))
    regex = r'_([0-9]_[0-9]+)'
    print('extract gene A id')
    dataframe['gene id: A'] = dataframe['gene A'].str.extract(regex, expand=True)
    print('extract gene B id')
    dataframe['gene id: B'] = dataframe['gene B'].str.extract(regex, expand=True)
    print('dataframe loading and node ID parsing time: {}'.format(datetime.now() - startTime))
    print('number of edges (dataframe rows): {}'.format(dataframe.shape[0]))
    return dataframe

def get_unique_nodes(df):
    df_nodes = df[['gene A', 'product A']].drop_duplicates()
    df_nodes.rename(columns={'gene A':'gene', 'product A':'product'}, inplace=True)
    df_nodes = \
        df_nodes.merge(df[['gene B', 'product B']].drop_duplicates().rename(
            columns={'gene B':'gene', 'product B':'product'}))
    return df_nodes

def add_nodes_from_df(network, df):
    # get set of unique nodes.
    startTime = datetime.now()
    df_nodes = get_unique_nodes(df)
    # load nodes
    print(df_nodes.head())
    for idx, row in df_nodes.iterrows():
        network.add_node(n=row['gene'], attr_dict={'product':row['product']})
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

def load_network(fname):
    df = load_df(fname)
    num_genes = len(set(df['gene A'].drop_duplicates().tolist() + df['gene B'].drop_duplicates().tolist()))
    print('number of unique nodes: {}'.format(num_genes))

    #n = networkx.DiGraph()  # directed
    n = networkx.Graph()   # undirected

    n = add_nodes_from_df(network=n, df=df)
    n = add_edges_from_df(network=n, df=df)

    return n
