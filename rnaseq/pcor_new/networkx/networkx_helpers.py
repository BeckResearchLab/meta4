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

def load_network_from_df(network, df):
    for idx, row in df.iterrows():
    # add_edge: The nodes u and v will be automatically added if
    # they are not already in the graph.
        network.add_edge(row['gene A'], row['gene B'],
                    attr_dict = {'pcor': row['pcor']} )
    return network

def load_network(fname):
    df = load_df(fname)
    num_genes = len(set(df['gene A'].drop_duplicates().tolist() + df['gene B'].drop_duplicates().tolist()))
    print('number of unique nodes: {}'.format(num_genes))

    startTime = datetime.now()
    #n = networkx.DiGraph()  # directed
    n = networkx.DiGraph()   # undirected
    network = load_network_from_df(n, df)
    print('networkx DiGraph filling time: {}'.format(datetime.now() - startTime))
    return network


