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
FILENAME = 'network.py.tsv'

def main():
    '''
    Extracts a co-occurence network from a pickle.

    Main entry point to code.
    '''

    # Read in the data
    print("reading previously saved data from pickle %s" % (DATA_PICKLE))
    with open(DATA_PICKLE, 'rb') as file:
        df = pickle.load(file)
        lwe = pickle.load(file)
        pmat = pickle.load(file)
        pcore_indices = pickle.load(file)
        pcor = pickle.load(file)
        lfdr_pcor = pickle.load(file)

    print("extracting off diagnol elements of precision matrix")
    pcor_indices = np.triu_indices(pmat.shape[0], 1)
    pcor = pmat[pcor_indices]

    # make dataframe network
    network = pd.DataFrame({"geneA":df.index.values[pcor_indices[0]], "geneB":df.index.values[pcor_indices[1]], "pcor":pcor})

    network.to_csv(FILENAME, sep='\t', index=False)

if __name__ == "__main__":
    main()
