import argparse
import glob
import os
import re
import sys

import pandas as pd


def prep_filename(gff_dir):
    # strip off trailing '/'
    gff_dir = gff_dir.rstrip('/')
    return os.path.join(gff_dir, '{}_gffs_concatenated.gff'.format(gff_dir))

def get_group_name(filepath):
    m = re.search(r'_group_([0-9]+)', filepath)
    if not m:
        print('failed to find gff batch number from {}'.format(filepath))
    return m.group(1)

def combine_gffs(gff_file_list, filename):
    """
    Combine a set of .gff files into a merged set.
    Do not include the coding sequences (encoded after ##FASTA)
    Make the locus names unique by replacing contigs_longer_than_XX_ to contigs_longer_than_XX_group_Z
    """
    # sort gff files, just in case
    gff_file_list = sorted(gff_file_list)

    with open(filename, 'w') as outfile:
        file_num = 1
        for f in gff_file_list:
            print('get the good stuff from {}'.format(f))
            group = get_group_name(f)
            with open(f) as f1:
                for line_num, line in enumerate(f1):        #keep the header from file1

                    # The first line is `##gff-version 3`.  Keep that for the first file.
                    if (file_num == 1) and (line_num == 0):
                        outfile.write(line)

                    # The end of the file has the entire FASTA sequence glued on.  Remove it.
                    elif '##FASTA' in line:
                        # We aren't keeping the nucleotide sequences!
                        break
                    # Delete subsequent lines like `##sequence-region k141_461591 1 2140`
                    elif (line_num > 0) and line.startswith('##'):
                        print('skip line: {}'.format(line))
                        continue
                    else:
                        # Need to give each file unique ID's. If not, each file has an ID like
                        # ID=contigs_longer_than_1500bp_00001
                        # and duplicates are a problem.
                        m = re.search('(contigs_[_a-z0-9]+bp)', line)
                        if m:
                            general_name = m.group(1)
                            new_name = general_name + '_group_' + group
                            line_edited = re.sub(general_name, new_name, line)
                            outfile.write(line_edited)
                        else:
                            print("failed to find re.search for '(contigs_[_a-z0-9]+bp)' in '{}'".format(line))
                            outfile.write(line)
                file_num += 1

def make_gene_table(gff_filename):
    colnames = ['contig', 'annotation tool', 'feature type', 'start', 'stop', '?', 'strand', '??', 'description']
    raw = pd.read_csv(gff_filename, sep='\t', names=colnames, skiprows=1)
    # E.g. D=contigs_longer_than_1500bp_group_1_00007;inference=ab initio ...
    raw['ID'] = raw['description'].str.extract('ID=([_A-z0-9]+);', expand=True)
    # Don't be as careful getting the gene products because we can go to the end of the line
    raw['product'] = raw['description'].str.extract('product=(.*)', expand=True)
    raw['bp'] = raw['stop'] - raw['start'] + 1
    # remove rows with contig == '##gff-version 3'
    raw = raw[ ~ raw['contig'].str.contains('##gff-version')]
    print(raw[raw['bp'].isnull()])
    raw['bp'] = raw['bp'].astype(int)
    raw[['contig', 'ID', 'product', 'bp']].to_csv(gff_filename + '.genes.tsv', sep='\t', index=False)



if __name__ == '__main__':
    # expect python 3 to be used.  Source activate py3
    assert sys.version_info > (3, 0), 'expected python 3; got {}'.format(sys.version)

    parser = argparse.ArgumentParser()
    parser.add_argument('gff_parent_folder', help='parent folder to get gffs from')
    args = parser.parse_args()

    if args.gff_parent_folder is None:
        print(args.print_help())
    else:
        print('gff parent folder: {}'.format(args.gff_parent_folder))

    gff_file_list = []
    search_path = args.gff_parent_folder + '/**/*.gff'
    print('look for all gff files at paths like {}'.format(search_path))
    for filename in glob.glob(search_path):
        gff_file_list.append(filename)
    print('found {} gff files: \n{}'.format(len(gff_file_list), gff_file_list))

    filename = prep_filename(args.gff_parent_folder)
    print('save concatenated gffs to {}'.format(filename))

    # combine the found gff files
    combine_gffs(gff_file_list, filename)
    make_gene_table(filename)

