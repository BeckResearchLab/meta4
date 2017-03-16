import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})

from plot_facets import bar_facets_from_pivoted_df


def plot_underscore_bars(input_df, filename=None, portrait=True):
    """
    """
    # make a copy, for safety.
    dataframe = input_df.copy()
    order_list = ['frac of RNA reads: __alignment_not_unique',
                  'frac of RNA reads: __ambiguous', 'frac of RNA reads: __no_feature',
                  'frac of RNA reads: __not_aligned', 'frac of RNA reads: __too_low_aQual']

    # specify colors, dark to light.
    colors = ['#810f7c', '#8856a7', '#8c96c6', '#9ebcda',    # purples
              '#006d2c', '#74c476'] # greens

    fig = plot_facets.bar_facets_from_pivoted_df(
        input_df=dataframe, pivoted=True,
        x='Genus italics', y='fraction of reads',
        order_list=order_list_italic, color_list=colors,
        portrait=portrait, filename=filename)
    fig = add_vline_to_all_subplots(fig, x=10.5-4, ymin=0, ymax=1, color='#636363')

    if filename is not None:
        fig.savefig(filename, bbox_inches='tight')

    return fig
