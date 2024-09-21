# this script is loading the fitting parameters from saved npy arrays and plot them versus age bins.

import numpy as np
import matplotlib
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import matplotlib.pyplot as plt
import sys

sim_type = np.array(['PE_PI_SN', 'noPE_PI_SN', 'PE_noPI_SN', 'PE_PI_SN_3Myr'])
age_bin = np.array(['10Myrs', '100Myrs', '5Myrs', '40Myrs'])


plt.clf()

for i in range (0, 4, 1):
        
        fit_single_10Myr = np.load('/home/nassim/analysis/clustering/fits/' + sim_type[i] + '/fit_single_10Myrs_' + sim_type[i] + '.npy')
        
        fit_broken_10Myr = np.load('/home/nassim/analysis/clustering/fits/' + sim_type[i] + '/fit_broken_10Myrs_' + sim_type[i] + '.npy')

        error_single_10Myr = np.load('/home/nassim/analysis/clustering/fits/' + sim_type[i] + '/error_single10Myrs_' + sim_type[i] + '.npy')
        
        error_broken_10Myr = np.load('/home/nassim/analysis/clustering/fits/' + sim_type[i] + '/error_broken10Myrs_' + sim_type[i] + '.npy')
        
        plt.figure(1)
        plt.subplot(121)
        plt.errorbar(np.arange(10, 100, 10), fit_broken_10Myr[:,2], yerr = error_broken_10Myr[:,2], label = sim_type[i])
        plt.subplot(122)
        plt.errorbar(np.arange(10, 100, 10), fit_broken_10Myr[:,0], yerr = error_broken_10Myr[:,0], label = sim_type[i])
        
        plt.figure(3)
        plt.errorbar(np.arange(10, 100, 10), fit_single_10Myr[:,0], yerr = error_single_10Myr[:,0], label = sim_type[i])
        

plt.figure(1)
plt.subplot(121)
plt.legend(loc="best")
plt.ylabel('Fitting slope')
plt.xlabel('Age bin')
plt.ylim(-2.0, 0)
plt.subplot(122)
plt.legend(loc="best")
plt.ylabel('Breaking point')
plt.xlabel('Age bin')
plt.ylim(0,4)
plt.savefig('/home/nassim/analysis/clustering/plots/compare/broken_fit_10Myr.png')

three = plt.figure(3)
plt.legend(loc="best")
plt.ylabel('Fitting slope')
plt.xlabel('Age bin')
plt.savefig('/home/nassim/analysis/clustering/plots/compare/single_fit_10Myr.png')



for i in range (0, 4, 1):
        
        fit_single_100Myr = np.load('/home/nassim/analysis/clustering/fits/' + sim_type[i] + '/fit_single_100Myrs_' + sim_type[i] + '.npy')
        
        fit_broken_100Myr = np.load('/home/nassim/analysis/clustering/fits/' + sim_type[i] + '/fit_broken_100Myrs_' + sim_type[i] + '.npy')

        error_single_100Myr = np.load('/home/nassim/analysis/clustering/fits/' + sim_type[i] + '/error_single100Myrs_' + sim_type[i] + '.npy')
        
        error_broken_100Myr = np.load('/home/nassim/analysis/clustering/fits/' + sim_type[i] + '/error_broken100Myrs_' + sim_type[i] + '.npy')
        
        plt.figure(2)
        plt.subplot(121)
        plt.errorbar(np.arange(100, 1000, 100), fit_broken_100Myr[:,2], yerr = error_broken_100Myr[:,2], label = sim_type[i])
        plt.subplot(122)
        plt.errorbar(np.arange(100, 1000, 100), fit_broken_100Myr[:,0], yerr = error_broken_100Myr[:,0], label = sim_type[i])

        plt.figure(4)
        plt.errorbar(np.arange(100, 1000, 100), fit_single_100Myr[:,0], yerr = error_single_100Myr[:,0], label = sim_type[i])
        
two = plt.figure(2, figsize = (20, 6))
plt.subplot(121)
plt.legend(loc="best")
plt.ylabel('Fitting slope')
plt.xlabel('Age bin')
plt.ylim(-1.5, 1.5)
plt.subplot(122)
plt.legend(loc="best")
plt.ylabel('Breaking point')
plt.xlabel('Age bin')
plt.ylim(0,4)
plt.savefig('/home/nassim/analysis/clustering/plots/compare/broken_fit_100Myr.png')

four = plt.figure(4)
plt.legend(loc="best")
plt.ylabel('Fitting slope')
plt.xlabel('Age bin')
plt.savefig('/home/nassim/analysis/clustering/plots/compare/single_fit_100Myr.png')
        
plt.show()