import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})

from collections import OrderedDict # for plot order

from plot_facets import bar_facets_from_pivoted_df
from plot_facets import add_vline_to_all_subplots


def plot_underscore_bars(input_df, filename=None, portrait=True):
    """
    """
    # make a copy, for safety.
    dataframe = input_df.copy()

    # prepare titles.
    # 'frac of RNA reads: __no_feature' --> 'no feature'
    dataframe.columns = [c.replace('frac of RNA reads: __', '') for c in dataframe.columns]
    dataframe.columns = [c.replace('_', ' ') for c in dataframe.columns]
    dataframe.rename(columns={'sum(frac RNA-seq mapped to genes)':'mapped to genes'}, inplace=True)

    colord = OrderedDict([
        ('mapped to genes','#31a354'), # dark green
        ('no feature', '#bdbdbd'), # light gray
        ('too low aQual', '#fa9fb5'), # pink
        ('not aligned', '#feb24c'), # light orange
        ('alignment not unique','w'),
        ('ambiguous','w'),
        ])

    fig = bar_facets_from_pivoted_df(
        input_df=dataframe, pre_pivoted=True, x='week',
        order_list=list(colord.keys()),
        color_list=list(colord.values()),
        portrait=portrait, filename=filename)
    fig = add_vline_to_all_subplots(fig, x=10.5-4, ymin=0, ymax=1, color='#636363')

    if filename is not None:
        fig.savefig(filename, bbox_inches='tight')

    return fig
