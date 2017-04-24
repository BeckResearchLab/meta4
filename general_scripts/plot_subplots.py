import matplotlib as mpl
# Seems to work better on AWS with Agg.
mpl.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

from functools import partial

# Need to use LaTeX to get italic fonts.
from matplotlib import rc
#rc('text', usetex=True)  # The italic O2 label stopped when I turned this off.
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
# http://stackoverflow.com/questions/2537868/sans-serif-math-with-latex-in-matplotlib
mpl.rcParams['text.latex.preamble'] = [
#       r'\usepackage{siunitx}',   # i need upright \micro symbols, but you need...
#       r'\sisetup{detect-all}',   # ...this to force siunitx to actually use your fonts
       r'\usepackage{helvet}',    # set the normal font here
       r'\usepackage{sansmath}',  # load up the sansmath so that math -> helvet
       r'\sansmath'               # <- tricky! -- gotta actually tell tex to use!
]

CONTIG_PARSE_REGEX = '[A-z0-9]+_([0-9]+_[0-9]+)'

REPLICATE_COLOR_DICT = {1:'#66c2a5', 2:'#fc8d62', 3:'#8da0cb', 4:'#e78ac3'}


def make_axd(axs, subplots, portrait=True):
    """
    axis dictionary for subplots.  Portrait or landscape orientation supported.
    """
    if subplots == 8:
        if portrait:
            return {('low', 1): axs[0, 0], ('high', 1): axs[0, 1],
                   ('low', 2): axs[1, 0], ('high', 2): axs[1, 1],
                   ('low', 3): axs[2, 0], ('high', 3): axs[2, 1],
                   ('low', 4): axs[3, 0], ('high', 4): axs[3, 1]}
        else:
             return {('low', 1): axs[0, 0], ('low', 2): axs[0, 1],
                     ('low', 3): axs[0, 2], ('low', 4): axs[0, 3],
                     ('high', 1): axs[1, 0], ('high', 2): axs[1, 1],
                     ('high', 3): axs[1, 2], ('high', 4): axs[1, 3]}
    elif subplots == 2:
        # handled the same, whether portrait or landscape.
        return {'low':axs[0], 'high':axs[1]}


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


def plot_multiple_series(df, groupby_var, color_dict, plot_function_partial, ax,
                         label_base=None):
    """
    Loop through all the series and plot.
    """
    for tup, plot_df in df.groupby(groupby_var):
        plot_df.sort_values('week')
        color = color_dict[tup]

        if label_base is not None:
            label = label_generator(label_base, tup)
        else:
            label = tup

        plot_function_partial(plot_df, color=color, ax=ax, label=label)


def plot_subplots(input_df, plot_function, ylabel, colors,
                  multiple_series=False, label_base=None, label=None,
                  portrait=True, subplots=8, legend_title=None,
                  prev_axs_and_axd=None):
    """
    General function to produce faceted plots (O2 vs rep or vice vers_dict)
    given already-pivotd data.

    :param input_df: a dataframe that may or may not be pivoted before plotting
    :param plot_function: plot function to use (likely a partial of a function)
    :param colors: list of hex (or RGB?) colors to use, or a dict of coolors if multiple_series=True
    :param portrait: True if tall is desired, or False if short
    :param subplots: 8 if one per oxygen/replicate combo, 2 if only separate by oxygen.
    :legend title: title to use for legend, if desired
    :return:
    """
    if portrait and prev_axs_and_axd is None:
        if subplots == 8:
            fig, axs = plt.subplots(4, 2, figsize=(10,10), sharex=True, sharey=True)
        elif subplots == 2:
            print('make it portrait')
            fig, axs = plt.subplots(2, 1, figsize=(3.5, 7), sharex=True, sharey=True)
            #plt.subplots_adjust(hspace=0.1)
        else:
            print('only 8 or 2 subplots are currently supported')
            return
    elif not portrait and prev_axs_and_axd is None:
        if subplots == 8:
            fig, axs = plt.subplots(2, 4, figsize=(14,8))
        elif subplots == 2:
            fig, axs = plt.subplots(1, 2, figsize=(7, 3.5))
            plt.subplots_adjust(wspace=0.5)
        else:
            print('only 8 or 2 subplots are currently supported')
            return
    if prev_axs_and_axd is None:
        axd = make_axd(axs, portrait=portrait, subplots=subplots)
    else:
        axs = prev_axs_and_axd[0]
        axd = prev_axs_and_axd[1]

    if subplots == 8 and portrait:
        axs[3,0].set_xlabel('week')
        axs[3,1].set_xlabel('week')
    elif subplots == 8 and not portrait:
        for x in [0, 1, 2, 3]:
                axs[0,x].set_xlabel('week')

    if subplots == 8:
        groupby_tuples = (['oxygen', 'replicate'])
    elif subplots == 2:
        groupby_tuples = ('oxygen')

    for tup, plot_df in input_df.groupby(groupby_tuples):
        oxygen_string = '$\mathregular{O_2}$'
        if subplots == 8:
            (o2, rep) = tup
            ax = axd[(o2, rep)]
            ax.set_title(o2 + ' ' + oxygen_string + ' replicate {}'.format(rep))
        elif subplots == 2:
            o2 = tup
            ax = axd[o2]
            #ax.set_title(o2 + ' $\mathregular{O_2}$')
            ax.set_title(o2 + ' ' + oxygen_string)
            ax.set_xlabel('week')

        #ax.set_title(o2 + ' oxygen' + ' replicate {}'.format(rep))
        plot_df.sort_values('week', inplace=True)

        if multiple_series:
            plot_multiple_series(plot_function_partial=plot_function, df=plot_df,
                                 ax=ax, groupby_var='replicate', color_dict=colors,
                                 label_base=label_base)
        else:
            # 170237 warning: adding label for scatter may break bars.
            if label is None:
                plot_function(plot_df=plot_df, ax=ax, colors=colors)
            else:
                plot_function(plot_df=plot_df, ax=ax, colors=colors, label=label)

    if portrait:
        # prevent subplot overlaps: set width, height to leave between subplots.
        plt.subplots_adjust(wspace = 0.1, hspace = 0.2)
        # add legend to the upper right
        if subplots == 8:
            axd[('high', 1)].legend(loc=(1.05, 0.4), title=legend_title)
        elif subplots == 2:
            axd['low'].legend(loc=(1.05, 0.2), title=legend_title)
    else:
        # prevent subplot overlaps: set width, height to leave between subplots.
        plt.subplots_adjust(wspace = 0.3, hspace = 0.3)
        # add legend to the upper right
        if subplots == 8:
            axd[('low', 4)].legend(loc=(1.05, 0.4), title=legend_title)
        elif subplots == 2:
            axd['high'].legend(loc=(1.05, 0.2), title=legend_title)

    # put labels up the left side.
    if subplots > 2:
        for ax in axs[:, 0]:
            ax.set_ylabel(ylabel)
    else:
        for ax in axs:
            ax.set_ylabel(ylabel)

    if prev_axs_and_axd is None:
        return fig
    else:
        return


def bar_subplots_from_pivoted_df(input_df, x, ylabel, order_list, color_list,
                               y=None, pre_pivoted=True, mulitple_series=False,
                               filename=None, portrait=True, legend_title=None):
    """

    :param order_list: order to plot columns by.  Should correspond to color list.
    :param pre_pivoted: True if dataframe was already pivoted with columns representing the desired bars
    """
    plot_function = partial(plot_bar, pre_pivoted=pre_pivoted, order_list=order_list)
    plot = plot_subplots(input_df, plot_function=plot_function, ylabel=ylabel,
                         colors=color_list, portrait=portrait,
                         legend_title=legend_title)
    return plot


def label_generator(base_string, string_to_replace_XX_with):
    return base_string.replace('XX', str(string_to_replace_XX_with))


def plot_scatter(plot_df, ax, x, y, colors, marker='o', linestyle='-', alpha=1,
                 label=None, label_generator_info=None):
    """
    :param colors: actually just one color in hex or RGB.
    # OLD:  a list of length 1.  This is to be consistent with the bar argument colors
    """
    color = colors
    ax.plot(plot_df[x], plot_df[y], color=color, marker=marker,
            linestyle=linestyle, label=label, alpha=alpha)


def plot_subplots_scatter(input_df, x, y, ylabel, marker, linestyle,
                          color_dict=REPLICATE_COLOR_DICT, portrait=True, subplots=2, label=''):
    plot_function = partial(plot_scatter, x=x, y=y, marker=marker, linestyle=linestyle)
    plot = plot_subplots(input_df, plot_function=plot_function, ylabel=ylabel, label=label,
                         colors=color_dict, portrait=portrait, multiple_series=True,
                         subplots=2, legend_title='replicate')
    return plot



