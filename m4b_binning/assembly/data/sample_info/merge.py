import pandas as pd

si = pd.read_csv('sample_info.tsv', sep='\t')
cr = pd.read_csv('meta4_sample_names--cryptic_to_sample_number.tsv', sep='\t')

m = pd.merge(cr, si, how='outer')
m.to_csv('sample_info_w_cryptic.tsv', sep='\t', index=False)
