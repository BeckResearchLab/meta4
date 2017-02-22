import argparse
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import re
import pandas as pd
import seaborn as sns

def read_tsv(f, cnames=None):
    if cnames is None:
        return pd.read_csv(f, sep='\t')
    else:
        return pd.read_csv(f, sep='\t', names=cnames)

def merge_two(one, two):
    m = pd.merge(one, two)
    assert m.shape[0] > 0, 'merge problem!  Got empty dataframe'
    return m

def extract_sample_name(df, col):
    regex = r'([0-9]+.[0-9]+.[0-9]+.[ACTG]+).*'
    return df[col].str.extract(regex, expand=True)

def merge_files(fq, bam, sample_info, cryptic_translator, type='metatranscriptome'):
    type_string = 'cryptic {} name'.format(type)

    fq = read_tsv(fq, ['fastq file', 'fastq reads']) 
    fq[type_string] = extract_sample_name(fq, 'fastq file') 

    bam = read_tsv(bam, ['bam file', 'bam reads'])
    bam[type_string] = extract_sample_name(bam, 'bam file') 

    si = read_tsv(sample_info)
    ct = read_tsv(cryptic_translator)
    # start merging. 
    results = merge_two(si, ct)
    results = merge_two(results, fq)
    results = merge_two(results, bam)
    return results

def plot_trend(dataframe):
    x = 'week'
    y1 = 'fastq reads'
    y2 = 'bam reads'
    
    fig, axs = plt.subplots(2, 4, figsize=(15, 6), sharex=True, sharey=True)
    axd = {('low', 1):axs[0, 0],
           ('low', 2):axs[0, 1],
           ('low', 3):axs[0, 2],
           ('low', 4):axs[0, 3], 
           ('high', 1):axs[1, 0],
           ('high', 2):axs[1, 1],
           ('high', 3):axs[1, 2],
           ('high', 4):axs[1, 3]}
    #colors = ['#66c2a5', '#fc8d62', '#8da0cb', '#e78ac3', '#a6d854']
    colors = ['#000000', '#31a354'] # black, green
    series = [y1, y2]
    colord = dict(zip(series, colors))

    for tup, df in dataframe.groupby(['oxygen', 'replicate']):
        ax = axd[tup]
        title = '{} O2, rep {}'.format(tup[0], tup[1])
        ax.set_title(title)
        df = df.copy()
        df.sort_values('week', ascending=False, inplace=True)
        for s in series:
            color = colord[s]
            ax.plot(df[x], df[s], color=color, label=s, marker='o', markersize=4)
        ax.set_xlabel(x)
        ax.set_ylabel('reads')
    ax.set_ylim(bottom=0)
            
    axs[0, 3].legend(bbox_to_anchor=(1.6, 1.))
    return fig


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('fastq_read_count_tsv', help='path to a tsv of read counts in a fastq file')
    parser.add_argument('bam_read_count_tsv', help='path to a tsv of read counts by bam file')
    parser.add_argument('merged_tsv_name', help='name of output file')
    args = parser.parse_args()
    print(args)

    if args.fastq_read_count_tsv is None:
        parser.print_help()
    else: 
        print('fastq_read_count_tsv: {}'.format(args.fastq_read_count_tsv))
        print('bam_read_count_tsv: {}'.format(args.bam_read_count_tsv))
        print('merged_tsv_name: {}'.format(args.merged_tsv_name))

    print('fastq_read_count_tsv: {}'.format(args.fastq_read_count_tsv))
    print('bam_read_count_tsv: {}'.format(args.bam_read_count_tsv))
    
    merged = merge_files(args.fastq_read_count_tsv, args.bam_read_count_tsv, 
                         '/work/m4b_binning/assembly/data/sample_info/sample_info.tsv', 
                         '/work/m4b_binning/assembly/data/sample_info/meta4_sample_names--cryptic_to_sample_number.tsv')

    merged.to_csv(args.merged_tsv_name, sep='\t', index=False)
    
    p = plot_trend(merged)
    p.savefig('summary--reads_in_fastq_and_bam.pdf', bbox_inches='tight')
    



