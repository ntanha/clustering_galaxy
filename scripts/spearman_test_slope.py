#this script is supposed to read the slope of the fits to the two point correlation function and calculate the spearman test

import numpy as np
import scipy
from scipy.stats import spearmanr

sim_type = np.array(['PE_PI_SN', 'noPE_PI_SN', 'PE_noPI_SN', 'PE_PI_SN_3Myr'])

for i in range (0, 4, 1):
        
        fit_single_100Myr = np.load('/home/nassim/analysis/clustering/fits/' + sim_type[i] + '/fit_single_100Myrs_' + sim_type[i] + '.npy')
        
        print 'spearmanr test for slopes of', sim_type[i], ' in 100Myrs age bins: ', scipy.stats.spearmanr(fit_single_100Myr[:,0], np.arange(100, 1000, 100))
