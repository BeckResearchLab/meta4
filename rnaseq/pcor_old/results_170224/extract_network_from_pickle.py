import os
import pickle
import sys

import matplotlib
matplotlib.use('Agg')

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

if os.path.isfile(DATA_PICKLE):
    print("reading previously saved data from pickle %s" % (DATA_PICKLE))
    with open(DATA_PICKLE, 'rb') as file:
        df = pickle.load(file)
        lwe = pickle.load(file)
        pmat = pickle.load(file)
        pcore_indices = pickle.load(file)
        pcor = pickle.load(file)
        lfdr_pcor = pickle.load(file)

