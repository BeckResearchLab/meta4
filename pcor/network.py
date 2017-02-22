#!/usr/bin/env python2.7
'''
Co-occurence network from expression data.
'''

import os
import pickle
import sys

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
import pandas as pd
import readline
from rpy2.robjects.packages import importr
from rpy2.robjects.vectors import FloatVector
from scipy import linalg
from sklearn.covariance import LedoitWolf

DATA_PICKLE = 'data.pkl'
FILENAME = 'normalized_counts.tsv'
PRUNE_GENES = 10000
PDF_FILENAME = 'network.py.pdf'

def main():
    '''
    Constructs a co-occurence network from gene expression data.

    Main entry point to code.
    '''

    # Read in the data
    if os.path.isfile(DATA_PICKLE):
        print("reading previously saved data from pickle %s" % (DATA_PICKLE))
        with open(DATA_PICKLE, 'rb') as file:
            df = pickle.load(file)
            lwe = pickle.load(file)
            pmat = pickle.load(file)
            pcore_indices = pickle.load(file)
            pcor = pickle.load(file)
            lfdr_pcor = pickle.load(file)
            #prob = pickle.load(file)
    else:
        print("reading in data from %s" % (FILENAME))
        df = pd.read_csv(FILENAME, sep='\t')
        print("found %d rows and %d columns" % (df.shape[0], df.shape[1]))
        # compute the row means and sort the data frame by descinding means
        df['row_means'] = df.mean(axis=1)
        df.sort_values('row_means', axis=0, ascending=False, inplace=True)
        df.drop('row_means', axis=1, inplace=True)
        # take the most abundant genes
        df = df.head(PRUNE_GENES)

        # Ledoit-Wolf optimal shrinkage coefficient estimate
        print("computing Ledoit-Wolf optimal shrinkage coeffecient estimate")
        lwe = LedoitWolf().fit(df.transpose())
        pmat = lwe.get_precision()
        # Convert symmetric matrix to array, first by getting indices
        # of the off diagonal elements, second by pulling them into
        # separate array (pcor).
        print("extracting off diagnol elements of precision matrix")
        pcor_indices = np.triu_indices(pmat.shape[0], 1)
        pcor = pmat[pcor_indices]

        # Determine edges by computing lfdr of pcor.
        print("computing lfdr of partial correlations")
        fdrtool = importr('fdrtool')
        lfdr_pcor = fdrtool.fdrtool(FloatVector(pcor), statistic="correlation", plot=False)
        #prob = 1-lfdr_pcor['lfdr']

        with open(DATA_PICKLE, 'wb') as file:
            pickle.dump(df, file, pickle.HIGHEST_PROTOCOL)
            pickle.dump(lwe, file, pickle.HIGHEST_PROTOCOL)
            pickle.dump(pmat, file, pickle.HIGHEST_PROTOCOL)
            pickle.dump(pcor_indices, file, pickle.HIGHEST_PROTOCOL)
            pickle.dump(pcor, file, pickle.HIGHEST_PROTOCOL)
            pickle.dump(lfdr_pcor, file, pickle.HIGHEST_PROTOCOL)
            #pickle.dump(prob, file, pickle.HIGHEST_PROTOCOL)

    print("making 1-lfdr vs. pcor plot")
    prob = 1-np.array(lfdr_pcor.rx2('lfdr'))
    with PdfPages(PDF_FILENAME) as pdf:
        plt.figure(figsize=(3, 3))
        plt.plot(range(7), [3, 1, 4, 1, 5, 9, 2], 'r-o')
        plt.title('Page One')
        pdf.savefig()  # saves the current figure into a pdf page
        plt.close()

        plt.plot(pcor[0:10000:10], prob[0:10000:10], 'o', markeredgecolor='k', markersize=3)
        plt.title("THIS IS A PLOT TITLE, YOU BET")
        plt.xlabel('partial correlation')
        plt.ylabel('lfdr')
        pdf.savefig
        plt.close()


if __name__ == "__main__":
    main()
