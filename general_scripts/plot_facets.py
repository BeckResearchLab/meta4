import matplotlib as mpl
# Seems to work better on AWS with Agg.
mpl.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

from functools import partial

# Need to use LaTeX to get italic fonts.
from matplotlib import rc
rc('text', usetex=True)
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
# http://stackoverflow.com/questions/2537868/sans-serif-math-with-latex-in-matplotlib
mpl.rcParams['text.latex.preamble'] = [
#       r'\usepackage{siunitx}',   # i need upright \micro symbols, but you need...
#       r'\sisetup{detect-all}',   # ...this to force siunitx to actually use your fonts
       r'\usepackage{helvet}',    # set the normal font here
       r'\usepackage{sansmath}',  # load up the sansmath so that math -> helvet
       r'\sansmath'               # <- tricky! -- gotta actually tell tex to use!
]


REPLICATE_COLOR_DICT = {1:'#66c2a5', 2:'#fc8d62', 3:'#8da0cb', 4:'#e78ac3'}


def axd_portrait(axs):
    """
    axis dictionary for portrait pages.
    """
    return {('low', 1): axs[0, 0], ('high', 1): axs[0, 1],
           ('low', 2): axs[1, 0], ('high', 2): axs[1, 1],
           ('low', 3): axs[2, 0], ('high', 3): axs[2, 1],
           ('low', 4): axs[3, 0], ('high', 4): axs[3, 1]}


def axd_landscape(axs):
    """
    axis dictionary for landscape pages.
    """
    return {('low', 1): axs[0, 0], ('low', 2): axs[0, 1], ('low', 3): axs[0, 2], ('low', 4): axs[0, 3],
            ('high', 1): axs[1, 0], ('high', 2): axs[1, 1], ('high', 3): axs[1, 2], ('high', 4): axs[1, 3]}


def add_vline_to_all_subplots(fig, x, ymin, ymax, color='#636363'):
    for ax in fig.axes:
        # It counts the x positions from 1, not from 4.
        ax.axvline(x=x, ymin=ymin, ymax=ymax, color=color)
    return fig


def plot_bar(plot_df, ax, pre_pivoted, colors, order_list):
    """
    :colors: a list of colors
    """

    if not pre_pivoted:
        assert y is not None, "need to specify fill value for pivot"
        plot_df = plot_df.pivot(index='week', columns=x, values=y)
    else:
        plot_df = plot_df.set_index('week')

    plot_df = plot_df[order_list]

    return plot_df.plot.bar(stacked=True, ax=ax, legend=False, color=colors)


def plot_multiple_series(df, groupby_var, color_dict, plot_function_partial, ax):
    """
    Loop through all the series and plot.
    """
    for tup, plot_df in df.groupby(groupby_var):
        plot_df.sort_values('week')
        color = color_dict[tup]
        plot_function_partial(plot_df, color=color, ax=ax)


def plot_facets(input_df, plot_function, colors, multiple_series=False, portrait=True, facets=8):
    """
    General function to produce faceted plots (O2 vs rep or vice vers_dict)
    given already-pivotd data.

    :param input_df: a dataframe that may or may not be pivoted before plotting
    :param plot_function: plot function to use
    :param colors: list of hex (or RGB?) colors to use, or a dict of coolors if multiple_series=True
    :param portrait: True if tall is desired, or False if short
    :return:
    """
    if portrait:
        fig, axs = plt.subplots(4, 2, figsize=(10,10))
        axd = axd_portrait(axs)
    else:
        fig, axs = plt.subplots(2, 4, figsize=(14,8))
        axd = axd_landscape(axs)

    for (o2, rep), plot_df in input_df.groupby(['oxygen', 'replicate']):
        ax = axd[(o2, rep)]
        ax.set_title(o2 + ' $\mathregular{O_2}$' + ' replicate {}'.format(rep))
        #ax.set_title(o2 + ' oxygen' + ' replicate {}'.format(rep))
        plot_df.sort_values('week', inplace=True)

        if multiple_series:
            plot_multiple_series(plot_function_partial=plot_function, df=plot_df,
                                 ax=ax, groupby_var='replicate', color_dict=colors)
        else:
            plot_function(plot_df=plot_df, ax=ax, colors=colors)

    if portrait:
        # prevent subplot overlaps: set width, height to leave between subplots.
        plt.subplots_adjust(wspace = 0.3, hspace = 0.6)
        # add legend to the upper right
        axd[('high', 1)].legend(loc=(1.05, 0))
    else:
        # prevent subplot overlaps: set width, height to leave between subplots.
        plt.subplots_adjust(wspace = 0.3, hspace = 0.3)
        # add legend to the upper right
        axd[('low', 4)].legend(loc=(1.05, 0.2))

    # put labels up the left side.
    for ax in axs[:, 0]:
        ax.set_ylabel('fractional abundance')

    return fig


def bar_facets_from_pivoted_df(input_df, x, order_list, color_list,
                               y=None, pre_pivoted=True, mulitple_series=False,
                               filename=None, portrait=True):
    """

    :param order_list: order to plot columns by.  Should correspond to color list.
    :param pre_pivoted: True if dataframe was already pivoted with columns representing the desired bars
    """
    plot_function = partial(plot_bar, pre_pivoted=pre_pivoted, order_list=order_list)
    plot = plot_facets(input_df, plot_function=plot_function, colors=color_list, portrait=True)
    return plot


def plot_scatter(plot_df, ax, x, y, color, marker='o', linestyle='-'):
    ax.plot(plot_df[x], plot_df[y], color=color, marker=marker, linestyle=linestyle)


def plot_facets_scatter(input_df, x, y, marker, linestyle, color_dict=REPLICATE_COLOR_DICT):
    plot_function = partial(plot_scatter, x=x, y=y, marker=marker, linestyle=linestyle)
    plot = plot_facets(input_df, plot_function=plot_function, colors=color_dict, portrait=True, multiple_series=True)
    return plot


