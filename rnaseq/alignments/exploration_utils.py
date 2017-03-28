import textwrap
import datetime

import matplotlib as mpl
import matplotlib.pyplot as plt

import itertools
import os
import re
import sys

import pandas as pd
import seaborn as sns

from functools import partial

sys.path.append('/work/general_scripts')
import plot_subplots  # 170327.  Most plots don't rely on this refactored code.

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
            '\(':'_',  # need to escape the (, ) chars
            '\)':'_',
            '%': '_'}
    for old, new in subs.items():
        name = re.sub(old, new, name)
    return name

def dict_of_sample_names():
    """
    Dict that converts from 8888.8.111111.GGAAGG type numbers to
    high_O2_replicate_3_week_2 type strings for PhD thesis.
    """
    si = si = pd.read_csv('/work/m4b_binning/assembly/data/sample_info/sample_info_w_cryptic.tsv',
                          sep='\t')
    si['name'] = si['oxygen'] + '_O2_rep_' + si['replicate'].astype(str) + "_week_" + si['week'].astype(str)
    #return si['cryptic metatranscriptome name'].tolist()
    return dict(zip(si['cryptic metatranscriptome name'].tolist(),
               si['name'].tolist()))


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

def load_frac_sums():
    counts = load_counts()
    frac_sums = counts.groupby('sample id')['frac RNA-seq reads'].sum()
    frac_sums_merged = pd.merge(et_sample_info(), frac_sums.to_frame().reset_index())
    return frac_sums

def load_counts_w_processing():
    """
    Takes a few minutes...
    """
    unders = load_underscore_stats()
    counts = load_counts()
    frac_sums = counts.groupby('sample id')['frac RNA-seq reads'].sum()
    frac_sums = frac_sums.reset_index().rename(
                        columns={'frac RNA-seq reads':'sum(frac RNA-seq mapped to genes)'})
    merged_df = pd.merge(unders, frac_sums, how='outer')

    # check sums
    merged_df['check sum'] = 0
    cnames_to_sum = [c for c in merged_df.columns
                     if (': __' in c) or (c == 'sum(frac RNA-seq mapped to genes)')]
    for c in cnames_to_sum:
        merged_df['check sum'] = merged_df['check sum'] + merged_df[c]

    return merged_df

def shorten_label(string, n):
    """
    Wrap lines of a long string, such that they are less than n characters wide.
    For axis labels.

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
        ax = axs_dict[o2]
        color = colors[rep]
        ax.plot(plot_df[x], plot_df[y],
                linestyle='-', marker='o', color=color, alpha = 0.8, label='rep {}'.format(rep))
    for ax_num, ax in enumerate(axs):
        ax.set_xlabel(x)
        if len(y) > 20:
            ylabel = shorten_label(y, 20)
            ax.set_ylabel(ylabel)
        else:
            ax.set_ylabel(y)
        facet_label = list(axs_dict.keys())[ax_num]
        ax.set_title('{} {}'.format(facet_label, facet_var))
    plt.tight_layout()
    axs[0].legend(bbox_to_anchor=(1.25, 1.))

    return fig


def prep_gene_cts(gene, counts, frac=True):
    """
    If frac=True, it's fraction relative to entire fastq counts, not to frac mapped.
    """
    sample_info = get_sample_info()
    if frac:
        value = 'frac RNA-seq reads'
    else:
        value = 'RNA reads'
    c = counts[counts['product'] == gene].groupby('sample id')[value].sum().reset_index()
    c = pd.merge(c, sample_info)
    return c

def plot_faceted_gene(df, colname):
    x = 'week'
    y = colname
    fig, axs = plt.subplots(2, 1, figsize=(6, 4), sharex=True, sharey=True)

    facet_var = 'oxygen'

    axs_dict = {'low': axs[0], 'high': axs[1]}
    colors = {1:'#66c2a5', 2:'#fc8d62', 3:'#8da0cb', 4:'#e78ac3'}

    for (facet_value, rep), plot_df in df.groupby([facet_var, 'replicate']):
        plot_df.sort_values('week', inplace=True)
        ax = axs_dict[facet_value]
        color = colors[rep]
        ax.plot(plot_df[x], plot_df[y],
                linestyle='-', marker='o', color=color, alpha = 0.8, label='rep {}'.format(rep))
        ax.set_title('{} {}'.format(facet_value, facet_var))
    for ax_num, ax in enumerate(axs):
        ax.set_xlabel(x)
        ax.set_ylabel('reads (RNA)')
    fig.suptitle(colname, fontsize=12)
    plt.subplots_adjust(top=0.85)
    axs[0].legend(bbox_to_anchor=(1.25, 1.))
    plt.subplots_adjust(hspace=0.45)

    return fig

def plot_read_counts_by_product(gene, sample_info):
    df = prep_gene_cts(gene, sample_info, frac=False)
    new_colname = 'RNA reads: {}'.format(gene)
    df.rename(columns={'RNA reads': new_colname}, inplace=True)
    p = plot_faceted_gene(df, new_colname)
    datestring = datetime.datetime.now().strftime("%y%m%d")
    fname = './figures/gene_reads/' + datestring + '_read_counts_' + filename_cleaner(gene) + '.pdf'
    print(fname)
    p.savefig(fname, bbox_inches='tight')

def plot_read_fracs_by_product(gene, sample_info, fignum=None):
    df = prep_gene_cts(gene, sample_info, frac=True)
    new_colname = 'frac RNA reads: {}'.format(gene)
    df.rename(columns={'frac RNA-seq reads': new_colname}, inplace=True)
    p = plot_faceted_gene(df, new_colname)
    folder = './figures/gene_read_fracs/'
    if not os.path.exists(folder):
        os.mkdir(folder)
    datestring = datetime.datetime.now().strftime("%y%m%d")
    if fignum is not None:
        fname = os.path.join(folder, datestring + '_read_fracs_' + str(fignum) + '_' + filename_cleaner(gene) + '.pdf')
    else:
        fname = os.path.join(folder, datestring + '_read_fracs_' + filename_cleaner(gene) + '.pdf')
    print(fname)
    p.savefig(fname, bbox_inches='tight')


def grab_genes_by_name(gene_name, dataframe):
    plot_df = dataframe[dataframe['product'] == gene_name].copy()
    plot_df['latex contig name'] = r'\mbox{' + plot_df['locus'] + r'}'
    plot_df['latex contig name'] = plot_df['latex contig name'].str.replace('_', '\_')
    return plot_df

def plot_abundance_of_genes_with_same_names(gene_name, dataframe, fignum=None, portrait=False, top_colors=None):
    plot_df = grab_genes_by_name(gene_name, dataframe)
    print("plot {} in each series' box".format(
        plot_df['locus'].drop_duplicates().shape[0]))
    x='week'
    y = 'frac RNA-seq reads'

    genes = plot_df.groupby('latex contig name')[y].max().sort_values(
        ascending=False).to_frame().reset_index()
    num_top=7
    top_genes = plot_df.groupby('latex contig name')[y].max().sort_values(
        ascending=False).to_frame().reset_index().head(num_top)['latex contig name'].tolist()

    if portrait:
        print('use portrait: 4 by 2')
        fig, axs = plt.subplots(4, 2, figsize=(10,10), sharex=True, sharey=True)
    else:
        fig, axs = plt.subplots(2, 4, figsize=(14,8), sharex=True, sharey=True)

    axd = plot_subplots.make_axd(axs, subplots=8, portrait=portrait)

    if top_colors is None:
        palette = itertools.cycle(sns.color_palette("Set1", num_top).as_hex())
    else:
        print('use specified colors: {}'.format(top_colors))
        palette = itertools.cycle(top_colors)

    top_gene_plot_fun =  partial(plot_subplots.plot_scatter, x=x, y=y,
                                 marker='o', linestyle='-')
    not_top_gene_plot_fun =  partial(plot_subplots.plot_scatter, x=x, y=y,
                                 marker='', linestyle='-')

    for locus, sub_df in plot_df.groupby('latex contig name'):
        if locus in top_genes:
            c = next(palette)
        sub_df = sub_df.copy()
        sub_df.sort_values('week', ascending=False, inplace=True)

        for tup, sub_sub_df in sub_df.groupby(['oxygen', 'replicate']):

            ax = axd[tup]
            if locus in top_genes:
                plot = plot_subplots.plot_subplots(sub_sub_df, plot_function=top_gene_plot_fun,
                            ylabel='fraction of fastq reads', colors=c, label=locus,
                            portrait=portrait, multiple_series=False, subplots=8, legend_title='',
                            prev_axs_and_axd = (axs, axd))
            else:
                c = '#909090'
                plot = plot_subplots.plot_subplots(sub_sub_df, plot_function=not_top_gene_plot_fun,
                            ylabel='fraction of fastq reads', label=locus, colors=c,
                            portrait=portrait, multiple_series=False, subplots=8, legend_title='',
                            prev_axs_and_axd = (axs, axd))
    fig.suptitle(gene_name, fontsize=16)
    plt.subplots_adjust(top=0.85)
    #for axv in [0, 1, 2, 3]:
    #    axs[1, axv].set_xlabel(x)

    folder = './figures/expression_by_locus/'
    if not os.path.exists(folder):
        os.mkdir(folder)
    datestring = datetime.datetime.now().strftime("%y%m%d")
    if portrait:
        orientation = 'portrait'
    else:
        orientation = 'landscape'
    if fignum is not None:
        fname = os.path.join(folder, datestring + '_loci_read_fracs_' + str(fignum) + '_' + filename_cleaner(gene_name) + '--' + orientation  + '.pdf')
    else:
        fname = os.path.join(folder, datestring + '_loci_read_fracs_' + filename_cleaner(gene_name) + '--' + orientation  +  '.pdf')

    fig.savefig(fname, bbox_inches='tight')
    return fig



