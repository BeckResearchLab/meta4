import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
# Need to use LaTeX to get italic fonts.
#rc('text', usetex=True)


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


def bar_facets_from_pivoted_df(input_df, x, order_list, color_list,
                               y=None, pre_pivoted=True,
                               filename=None, portrait=True):
    """
    General function to produce faceted plots (O2 vs rep or vice versa)
    given already-pivotd data.

    :param input_df: a dataframe that may or may not be pivoted before plotting
    :param pre_pivoted: True if dataframe was already pivoted with columns representing the desired bars 
    :param x: value to use as x in pivot (and plot)
    :param y: value to use as y (bar rectangle) in pivot and plot
    :param order_list: order to plot columns by.  Should correspond to color list.
    :param color_list: list of hex (or RGB?) colors to use
    :param filename: filename to save the file with
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
        plot_df.sort_values('week', inplace=True)

        if not pre_pivoted:
            assert y is not None, "need to specify fill value for pivot"
            plot_df = plot_df.pivot(index='week', columns=x, values=y)
        else:
            plot_df = plot_df.set_index('week')
        
        plot_df = plot_df[order_list]
        plot_df.plot.bar(stacked=True, ax=ax, legend=False, color=color_list)

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



























