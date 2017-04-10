import argparse
import sys

sys.path.append('../')
sys.path.append('/work/rnaseq/pcor_new')

from network_prep import load_and_drop_rare_features

def load_and_drop_rare(min_percent, fname):
    df = load_and_drop_rare_features(min_percent=min_percent, plot_dist=False)
    df.to_csv(fname, sep='\t')

def load_without_drop(fname):
    # NEW: Don't drop rare features.
    # Actually kind of wasteful
    print('load data, in format percent of fastq')
    df = load_and_drop_rare_features(plot_dist=False)
    print('write un-trimmed out put to {}'.format(fname))
    df.to_csv(fname, sep='\t')

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("--percent", type=float,
                        help="feature reduction parameter: percent of reads in the fastq for at least one sample")
    parser.add_argument("--dest_filename", type=str,
                        help="filename to save trimmed files at")
    args = parser.parse_args()
    print('arguments have been parsed')

    # check that the arguments were supplied.
    try:
        options = parser.parse_args()
        assert args.dest_filename is not None, \
            print('need to specify the destination filename')
    except:
        parser.print_help()
        sys.exit(0)

    print('save input to R as ', args.dest_filename)

    if args.percent:
        load_and_drop_rare(args.percent, args.dest_filename)
    else:
        load_without_drop(args.dest_filename)

