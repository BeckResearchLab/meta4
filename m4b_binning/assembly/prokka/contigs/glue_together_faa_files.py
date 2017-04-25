import argparse
import glob
import os
import re
import sys

import pandas as pd


def prep_filename(faa_dir):
    # strip off trailing '/'
    faa_dir = faa_dir.rstrip('/')
    return os.path.join(faa_dir, '{}_faa_files_concatenated.faa'.format(faa_dir))

def get_group_name(filepath):
    m = re.search(r'_group_([0-9]+)', filepath)
    if not m:
        print('failed to find faa batch number from {}'.format(filepath))
    return m.group(1)

def combine_faa_files(faa_file_list, filename):
    """
    Just glue them together, but change strings like contigs_longer_than_XX_ to contigs_longer_than_XX_group_Z
    Matches gff glue strategy.
    """
    # sort faa files, just in case
    faa_file_list = sorted(faa_file_list)

    with open(filename, 'w') as outfile:
        for f in faa_file_list:
            print('get the good stuff from {}'.format(f))
            group = get_group_name(f)
            with open(f) as f1:
                for line_num, line in enumerate(f1):        #keep the header from file1
                    # Need to give each gene a unique ID that matches the .gff file.
                    # The original strings are like contigs_longer_than_1500bp_00001
                    # and duplicates across Prokka calls are a problem.
                    m = re.search('(contigs_[_a-z0-9]+bp)', line)
                    if m:
                        general_name = m.group(1)
                        new_name = general_name + '_group_' + group
                        line_edited = re.sub(general_name, new_name, line)
                        outfile.write(line_edited)
                    else:
                        outfile.write(line)

if __name__ == '__main__':
    # expect python 3 to be used.  Source activate py3
    assert sys.version_info > (3, 0), 'expected python 3; got {}'.format(sys.version)

    parser = argparse.ArgumentParser()
    parser.add_argument('faa_parent_folder', help='parent folder to get faa files from')
    args = parser.parse_args()

    if args.faa_parent_folder is None:
        print(args.print_help())
    else:
        print('faa parent folder: {}'.format(args.faa_parent_folder))

    faa_file_list = []
    search_path = args.faa_parent_folder + '/**/*.faa'
    print('look for all faa files at paths like {}'.format(search_path))
    for filename in glob.glob(search_path):
        faa_file_list.append(filename)
    print('found {} faa files: \n{}'.format(len(faa_file_list), faa_file_list))

    filename = prep_filename(args.faa_parent_folder)
    print('save concatenated faa files to {}'.format(filename))

    # combine the found faa files
    combine_faa_files(faa_file_list, filename)

