import textwrap

import matplotlib as mpl
import matplotlib.pyplot as plt

import re

import pandas as pd
import seaborn as sns

def load_cryptic_to_sample_names():
    return pd.read_csv('/work/m4b_binning/assembly/data/sample_info/meta4_sample_names--cryptic_to_sample_number.tsv', sep='\t')

def load_rna_fastq_read_counts():
    counts = pd.read_csv('/work/rnaseq/rna_fastq_read_counts.tsv', sep='\t', names=['fastq', 'RNA reads in fastq'])
    # trim off fastq.gz
    counts['cryptic metatranscriptome name'] =  counts['fastq'].str.extract(r'(.*).fastq.gz')
    csn = load_cryptic_to_sample_names()
    si = get_sample_info()
    df = pd.merge(csn, si, how='outer')
    df = pd.merge(df, counts)
    return df

def load_unders():
    # Load the underscore portion of the data. 
    unders = pd.read_csv(
        #'./map_to_contigs_longer_than_1500bp/contigs_longer_than_1500bp_melted_underscores.tsv', 
        './map_to_contigs_longer_than_1500bp/map_to_contigs_longer_than_1500bp_melted_underscores.tsv', 
        sep='\t')
    return unders

def get_sample_info():
    unders = load_unders()
    sample_info = unders[['sample id', 'oxygen', 'replicate', 'week']].drop_duplicates()
    return sample_info

def filename_cleaner(name):
    subs = {'/': '--',
            ' ': '_',
            ':': '',
            '%': '_'}
    for old, new in subs.items():
        name = re.sub(old, new, name)
    return name


def load_underscore_stats(rna_fastq=True):
    unders = load_unders()
    unders_pivoted = unders.pivot(index='sample id', columns='locus', values='RNA reads').reset_index()
    unders_pivoted = pd.merge(unders_pivoted, get_sample_info())
    if rna_fastq:
        print('merge column informing total reads in fastq')
        fq = load_rna_fastq_read_counts()
        unders_pivoted = pd.merge(unders_pivoted, fq, how='left')
        if 'fastq' in unders_pivoted.columns:
            del unders_pivoted['fastq']
        underscore_cols = [c for c in unders_pivoted.columns if '__' in c]
        print('underscore_cols: {}'.format(underscore_cols))
        for cname in underscore_cols:
            c_frac_name = 'frac of RNA reads: {}'.format(cname)
            unders_pivoted[c_frac_name] = unders_pivoted[cname]/unders_pivoted['RNA reads in fastq']
        
    return unders_pivoted

def load_counts(rna_fastq=True):
    #counts = pd.read_csv('./map_to_contigs_longer_than_1500bp/contigs_longer_than_1500bp_melted.tsv', sep='\t')
    counts = pd.read_csv('./map_to_contigs_longer_than_1500bp/map_to_contigs_longer_than_1500bp_melted.tsv', sep='\t')
    if rna_fastq:
        fq = load_rna_fastq_read_counts()
        counts = pd.merge(counts, fq)
        # Make a column for the fraction of RNA-seq reads each gene reperesents.
        # It is possible these do not sum to 1. 
        counts['frac RNA-seq reads'] = counts['RNA reads']/counts['RNA reads in fastq']
    return counts

def shorten_label(string, n):
    """
    E.g. 'frac of RNA reads: __no_feature', 20 --> 2 lines:
        frac of RNA reads:
        __no_feature
    """
    return textwrap.fill(string, n)
   
def plot_faceted(df, colname):
    x = 'week'
    y = colname
    fig, axs = plt.subplots(2, 1, figsize=(6, 4), sharex=True, sharey=True)
    
    facet_var = 'oxygen'
    
    axs_dict = {'low': axs[0], 'high': axs[1]}
    colors = {1:'#66c2a5', 2:'#fc8d62', 3:'#8da0cb', 4:'#e78ac3'}
    
    for (o2, rep), plot_df in df.groupby([facet_var, 'replicate']):
        plot_df.sort_values('week', inplace=True)
        #print(df.head())
        ax = axs_dict[o2]
        #print(ax)
        color = colors[rep]
        #print(plot_df[y].head())
        ax.plot(plot_df[x], plot_df[y],
                linestyle='-', marker='o', color=color, alpha = 0.8, label='rep {}'.format(rep))
        #ax.plot(plot_df['week'], plot_df['sum'], linestyle='-', marker='o', color=color)
    for ax_num, ax in enumerate(axs):
        ax.set_xlabel(x)
        if len(y) > 20:
            ylabel = shorten_label(y, 20) 
            ax.set_ylabel(ylabel)
        else:
            ax.set_ylabel(y)
        #print(list(axs_dict.keys()))
        facet_label = list(axs_dict.keys())[ax_num]
        ax.set_title('{} {}'.format(facet_label, facet_var))
    plt.tight_layout()
    axs[0].legend(bbox_to_anchor=(1.25, 1.))
    
    return fig


def prep_gene_cts(gene, counts):
    """
    If frac=True, it's fraction relative to entire fastq counts, not to frac mapped. 
    """
    sample_info = get_sample_info()
    if frac:
        value = 'frac RNA-seq reads'
    else:
        value = 'RNA reads'
    c = counts[counts['product'] == gene].groupby('sample id')[value].sum().reset_index()
    print(c.shape)
    print(c.head(2))
    print(c.columns)
    c = pd.merge(c, sample_info)
    return c
    #sample_info = get_sample_info()
    #c = counts[counts['product'] == gene].groupby('sample id')['RNA reads'].sum().reset_index()
    #print(c.shape)
    #print(c.head(2))
    #print(c.columns)
    #c = pd.merge(c, sample_info)
    #return c

def plot_faceted_gene(df, colname):
    x = 'week'
    y = colname
    fig, axs = plt.subplots(2, 1, figsize=(6, 4), sharex=True, sharey=True)
    
    facet_var = 'oxygen'
    
    axs_dict = {'low': axs[0], 'high': axs[1]}
    colors = {1:'#66c2a5', 2:'#fc8d62', 3:'#8da0cb', 4:'#e78ac3'}
    
    for (facet_value, rep), plot_df in df.groupby([facet_var, 'replicate']):
        plot_df.sort_values('week', inplace=True)
        #print(df.head())
        ax = axs_dict[facet_value]
        #print(ax)
        color = colors[rep]
        #print(plot_df[y].head())
        ax.plot(plot_df[x], plot_df[y],
                linestyle='-', marker='o', color=color, alpha = 0.8, label='rep {}'.format(rep))
        #ax.plot(plot_df['week'], plot_df['sum'], linestyle='-', marker='o', color=color)
        ax.set_title('{} {}'.format(facet_value, facet_var))
    for ax_num, ax in enumerate(axs):
        ax.set_xlabel(x)
        ax.set_ylabel('reads (RNA)')
    fig.suptitle(colname, fontsize=12)
    plt.subplots_adjust(top=0.85)
    #plt.tight_layout()
    axs[0].legend(bbox_to_anchor=(1.25, 1.))
    plt.subplots_adjust(hspace=0.45)
    
    return fig

def plot_read_counts_by_product(gene, sample_info):
    df = prep_gene_cts(gene, sample_info)
    new_colname = 'RNA reads: {}'.format(gene)
    df.rename(columns={'RNA reads': new_colname}, inplace=True)
    p = plot_faceted_gene(df, new_colname)
    fname = './figures/gene_reads/170213_read_counts_' + filename_cleaner(gene) + '.pdf'
    print(fname)
    p.savefig(fname, bbox_inches='tight')
