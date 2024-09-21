#this script reads into npy files produced by mst_hist.py script and make histograms with logarythmic scale

import numpy as np
import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import cm

low = np.log10(2.)
high = np.log10(1000.)
bins = 50
width = (high - low)/bins

bins_j = np.load('bins.npy')
average_hist_less40 = np.load('average_hist_less40.npy')

plt.bar(bins_j[0:-1], average_hist_less40, width,  align='edge', log=True)
axes = plt.gca()
axes.set_ylim([0,60])
plt.show()