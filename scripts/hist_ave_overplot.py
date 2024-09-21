# this script loads the average_hist_less40_90t0990.npy produced by the mst_hist.py script from the different simulations and overplots them. 

import numpy as np
import matplotlib
import matplotlib.pyplot as plt

"""loading the data"""
bins = np.load("/home/nassim/analysis/clustering/histograms/noPE_PI_SN/bins.npy")                                                  #position of the logarythmic bins
PE_PI_SN = np.load("/home/nassim/analysis/clustering/histograms/PE_PI_SN/average_hist_less40_90to990.npy")                            #average histogram of the simulation withough photo-electric heating
noPE_PI_SN = np.load("/home/nassim/analysis/clustering/histograms/noPE_PI_SN/average_hist_less40_90to990.npy")                        #average histogram of the simulation with all of the  processes
PE_noPI_SN = np.load('/home/nassim/analysis/clustering/histograms/PE_noPI_SN/average_hist_less40_90to990.npy')
PE_PI_SN_3Myr = np.load('/home/nassim/analysis/clustering/histograms/PE_PI_SN_3Myr/average_hist_less40_90to990.npy')

bin_num = np.size(bins)
low = np.log10(2.)
high = np.log10(1000.)
width = (high - low)/bin_num

"""plotting"""
plt.clf()
plt.figure(1)
plt.plot(bins[0:-1], PE_PI_SN, markersize=2, label='PE_PI_SN')
plt.plot(bins[0:-1], noPE_PI_SN, markersize=2, label='noPE_PI_SN')
plt.plot(bins[0:-1], PE_noPI_SN, markersize=2, label='PE_noPI_SN')
plt.plot(bins[0:-1], PE_PI_SN_3Myr, markersize=2, label='PE_PI_SN_3Myr')
plt.xlabel('seperation [log pc]', fontsize=15)
plt.ylabel('normalized frequency', fontsize=15)
axes = plt.gca()
axes.set_ylim([0,1])
plt.legend(loc="best")
plt.savefig("../plots/compare/overplot_hist_average_less40_90to990.png")

plt.figure(2)
plt.bar(bins[0:-1], PE_PI_SN, width,  align='edge', color = 'white', alpha=1, edgecolor = 'blue', label='PE_PI_SN')
plt.bar(bins[0:-1], noPE_PI_SN, width,  align='edge', color = 'white', alpha=0.5, edgecolor = 'orange', label='noPE_PI_SN')
plt.bar(bins[0:-1], PE_noPI_SN, width,  align='edge', color = 'white', alpha=0.5, edgecolor = 'green', label='PE_noPI_SN')
plt.bar(bins[0:-1], PE_PI_SN_3Myr, width,  align='edge', color = 'white', alpha=0.5, edgecolor = 'red', label='PE_PI_SN_3Myr')
plt.savefig('../plots/compare/overplot_hist_average_less40_90to990_box.png')
plt.show()
