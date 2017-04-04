import pandas as pd

def make_merged_df():
    si = pd.read_csv('sample_info.tsv', sep='\t')
    cr = pd.read_csv('meta4_sample_names--cryptic_to_sample_number.tsv', sep='\t')

    m = pd.merge(cr, si, how='outer')
    return m

def prep_sample_desc_string(merged_df, tex=False):
    if tex:
        O2 = ' $\mathregular{O_2}$ '
    else:
        O2 = ' O2 '
    merged_df['descriptive string'] = merged_df['oxygen'] + O2 + 'week ' +  merged_df['week'].astype(str) + ' rep ' +  merged_df['replicate'].astype(str)
    return merged_df


if __name__ == '__main__':
    m = make_merged_df()
    m = prep_sample_desc_string(m, tex=True)
    print(m.head(3))
    m.to_csv('sample_info_w_cryptic.tsv', sep='\t', index=False)


