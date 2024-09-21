# this is a script to overplot the mean, sigma and the skewness of the mst for different simulation types

import numpy as np

import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm

sim_type = np.array(['PE_PI_SN', 'noPE_PI_SN', 'PE_noPI_SN', 'PE_PI_SN_3Myr'])

t = np.arange(90, 990, 10)
plt.clf()
for i in range (0, 4, 1):
    l_conv = np.load('/home/nassim/analysis/clustering/mst_stat/' + sim_type[i] + '/l_conv.npy')
    sigma_conv = np.load('/home/nassim/analysis/clustering/mst_stat/' + sim_type[i] + '/sigma_conv.npy')
    gamma_conv = np.load('/home/nassim/analysis/clustering/mst_stat/' + sim_type[i] + '/gamma_conv.npy')
    sfr_conv = np.load('/home/nassim/analysis/clustering/mst_stat/' + sim_type[i] + '/sfr_conv.npy')
    
    plt.figure(1)
    plt.plot(t, l_conv, markersize=1.5, label = sim_type[i])

    plt.figure(2)
    plt.plot(t, sigma_conv, markersize=1.5, label = sim_type[i])

    plt.figure(3)
    plt.plot(t, gamma_conv, markersize=1.5, label = sim_type[i])
    
    plt.figure(4)
    plt.plot(t, sfr_conv, markersize=1.5, label = sim_type[i])
    

plt.figure(1)
plt.ylabel('length mean[pc]', fontsize=12)
plt.xlabel('time[Myr]', fontsize=12)
plt.legend(loc = 'best')
plt.savefig('../plots/compare/length_mean.png')      
            
plt.figure(2)
plt.ylabel('sigma[pc]', fontsize=12)
plt.xlabel('time[Myr]', fontsize=12)
plt.legend(loc = 'best')
plt.savefig('../plots/compare/length_sigma.png')

plt.figure(3)
plt.ylabel('skewness', fontsize=12)
plt.xlabel('time[Myr]', fontsize=12)
plt.legend(loc = 'best')
plt.savefig('../plots/compare/length_skewness.png')

plt.figure(4)
plt.ylabel('SFR', fontsize=12)
plt.xlabel('time[Myr]', fontsize=12)
plt.legend(loc = 'best')
plt.savefig('../plots/compare/sfr.png')




        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        