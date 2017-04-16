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

def draw(G, node_label='product', layout=nx.circular_layout, edge_multiplier=100):
    node_labels = nx.get_node_attributes(G, node_label) # usually label by gene product

    edge_labels = nx.get_edge_attributes(G, 'pcor')
    edge_labels = round_values_in_dict(edge_labels, n=2)
    edge_widths = [ d['pcor']*edge_multiplier for (u,v,d) in G.edges(data=True)]
    def assign_color(pcor):
        if pcor > 0:
            return "#1b9e77" # green
        else:
            return "#d95f02" # orange
    edge_colors = [assign_color(d['pcor']) for (u,v,d) in G.edges(data=True)]

    pos = layout(G)
    nx.draw(G, pos, width=edge_widths, edge_color=edge_colors)
    nx.draw_networkx_labels(G, pos,
                            labels=node_labels,
                            font_size=9)

    nx.draw_networkx_edge_labels(G, pos,
                                 edge_labels=edge_labels,
                                 font_size=9)


def subgraph_by_cutoff(G, cutoff, hypothetical=True):
    # trim out hypotheticals *first*
    if not hypothetical:
        non_hypo_tuples = [n for n in G.nodes_iter(data=True)
                           if 'hypothetical' not in n[1]['product']]
        non_hypo_nodes = [ID for (ID, node_dict) in non_hypo_tuples]
        SG = nx.Graph(G.subgraph(non_hypo_nodes))
    else:
        SG = G.copy()

    # Next, make the sub-graph.  This loses node attributes, so we upgrae it below.
    SG0 = nx.Graph([(u,v,d) for u,v,d in SG.edges(data=True)  # get nodes, and edge attributes.
                   if abs(d['pcor']) > cutoff], strict=True)
    # upgrade to a graph with node attributes
    SG = nx.Graph(SG.subgraph(SG0.nodes()))

    # trim out hypotheticals

    return SG

def get_nodes_including_string(G, string):
    return [n for n in G.nodes_iter(data=True)
            if string in n[1]['product']]

