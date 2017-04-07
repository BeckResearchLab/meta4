import sys

sys.path.append('../')

from network_prep import load_and_drop_rare_features

MIN_PERCENT=0.005

df = load_and_drop_rare_features(min_percent=0.005, plot_dist=False)
df.to_csv('input_for_R--min_percent_{}--unnormalized.tsv'.format(MIN_PERCENT), sep='\t')


