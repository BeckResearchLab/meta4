import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})

from plot_facets import bar_facets_from_pivoted_df
from plot_facets import add_vline_to_all_subplots


def plot_underscore_bars(input_df, filename=None, portrait=True):
    """
    """
    # make a copy, for safety.
    dataframe = input_df.copy()
    order_list = ['sum(frac RNA-seq mapped to genes)',
                  'frac of RNA reads: __alignment_not_unique',
                  'frac of RNA reads: __ambiguous', 'frac of RNA reads: __no_feature',
                  'frac of RNA reads: __not_aligned', 'frac of RNA reads: __too_low_aQual']

    # specify colors
    colors = ['#31a354', # dark green,
              'w',
              'w',
              '#bdbdbd', # light gray
              #'#f03b20', # dark orange
              '#feb24c', # light orange
              '#fa9fb5', # pink
             ]

    fig = bar_facets_from_pivoted_df(
        input_df=dataframe, pre_pivoted=True, x='week',
        order_list=order_list, color_list=colors,
        portrait=portrait, filename=filename)
    fig = add_vline_to_all_subplots(fig, x=10.5-4, ymin=0, ymax=1, color='#636363')

    if filename is not None:
        fig.savefig(filename, bbox_inches='tight')

    return fig
