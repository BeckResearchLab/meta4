import matplotlib as mpl
import matplotlib.pyplot as plt
import itertools  # for color palette cycling
import re
import pandas as pd
import seaborn as sns
import sys
from cycler import cycler
import seaborn as sns

import networkx as nx
if sys.version_info.major == 2:
    import community
else:
    print("cannot use community detection algorithms in pyton 3 (that's you!)")

def plot_degree_rank(G):
    degree_sequence=sorted(nx.degree(G).values(),reverse=True) # degree sequence
    #print "Degree sequence", degree_sequence
    dmax=max(degree_sequence)

    fig, ax = plt.subplots(1,1, figsize=(4, 2.5))
    ax.loglog(degree_sequence,'b-',marker='o', color='#43a2ca')
    plt.title("Degree rank plot")
    plt.ylabel("degree")
    plt.xlabel("rank")
    return fig

def graph_by_nodes_list(G, nodes_list):
    return G.subgraph(nodes_list)

def round_values_in_dict(d, n=2):
    for k, v in d.items():
        d[k] = round(v, n)
    return d

def draw(G, node_label='product', layout=nx.circular_layout):
    node_labels = nx.get_node_attributes(G, node_label) # usually label by gene product

    edge_labels = nx.get_edge_attributes(G, 'pcor')
    edge_labels = round_values_in_dict(edge_labels, n=2)

    pos = layout(G)
    nx.draw(G, pos)
    nx.draw_networkx_labels(G, pos,
                                  labels=node_labels,
                                  font_size=9)

    nx.draw_networkx_edge_labels(G, pos,
                                       edge_labels=edge_labels,
                                       font_size=9)


