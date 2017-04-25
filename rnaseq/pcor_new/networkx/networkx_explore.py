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

def draw(G, node_label='product', layout=nx.circular_layout,
         edge_multiplier=100):
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
    nx.draw(G, pos, node_size=15, width=edge_widths, edge_color=edge_colors)

    print('include node labels')
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
    tuple_list =  [n for n in G.nodes_iter(data=True)
                   if string in n[1]['product']]
    return tuple_list

def get_nodes_including_list_of_strings(G, string_list):
    nodes_tuples = []
    for s in string_list:
        nodes_tuples = nodes_tuples + get_nodes_including_string(G, s)
        #lists.append(get_nodes_including_string(G, s))
    node_IDs = [ID for (ID, node_dict) in nodes_tuples]
    unique_node_IDs = list(set(node_IDs))
    print(unique_node_IDs)
    print('of {} nodes, {} were unique'.format(len(node_IDs), len(unique_node_IDs)))
    return unique_node_IDs


def sub_abs_pocr_df(G):
    def sum_abs_pcors_by_node(G, node):
        return sum([abs(d['pcor']) for n, d in G[node].items()])

    sub_abs_dict = {n: sum_abs_pcors_by_node(G, n) for n, d in G.degree().items()}
    node_info = pd.DataFrame.from_dict(sub_abs_dict, orient='index')
    node_info.rename(columns={0: "sum(abs(pcor))"}, inplace=True)
    node_info.sort_values('sum(abs(pcor))', ascending=False, inplace=True)
    return node_info

def node_names_df(G):
    n_dict = {n:d['product'] for (n, d) in G.nodes(data=True)}
    n_df = pd.DataFrame.from_dict(n_dict, orient='index')
    n_df.rename(columns={0: "product"}, inplace=True)
    return n_df

def num_edges_df(G):
    n_edges_dict = {n:len(G[n]) for n in G.nodes()}
    node_info = pd.DataFrame.from_dict(n_edges_dict, orient='index')
    node_info.rename(columns={0: "# edges"}, inplace=True)
    node_info.sort_values("# edges", ascending=False, inplace=True)
    return node_info

def summarise_sum_abs_pcors(G):
    pcors = sub_abs_pocr_df(G).reset_index()
    names = node_names_df(G).reset_index()
    merged = pd.merge(names, pcors, how='outer')
    merged.sort_values('sum(abs(pcor))', ascending=False,
                       inplace=True)

    num_edges = num_edges_df(G).reset_index()
    merged = pd.merge(merged, num_edges, how='outer')
    return merged


