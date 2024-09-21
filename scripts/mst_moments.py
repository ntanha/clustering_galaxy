# this script reads into dd-array ascii files, calculate the minimum spanning tree and calculated the first, second and third moment as m_k = (l - l_mean)^k

import numpy as np
import scipy
from scipy.optimize import curve_fit
from scipy.sparse.csgraph import minimum_spanning_tree
from scipy.stats import skewtest
from scipy.stats import spearmanr
import math
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm


jump = 1
l_mean = np.zeros((990-90)/10)
moment_1 = np.zeros((990-90)/10)
moment_2 = np.zeros((990-90)/10)
moment_3 = np.zeros((990-90)/10)
sigma = np.zeros((990-90)/10)							#standard deviation
gamma = np.zeros((990-90)/10)							#skewness

i = 0
for k in range (90, 990, 10):
    kk = k*jump
    if (kk < 10):
	num = '00' + str(kk)
    elif (kk < 100):
	num = '0' + str(kk)
    elif (kk < 1000):
	num = str(kk)
	
    filename_path = "/home/nassim/analysis/clustering/dd_less40/pe_pi_sn"
    filename_base = filename_path + "/dd_less40_"
    filename = filename_base + num + ".txt"
    
    print ''
    print 'Read dd-array file ', num
    
    FullTree = np.genfromtxt(filename, dtype=float)
    mst = minimum_spanning_tree(FullTree)
    
    
    Nstar = np.shape(mst)[0]
    print "number of stars in this snapshot between age 0 and 40 Myrs is %d" %Nstar
    
    L = np.zeros(Nstar-1)							#an array of all the non-zero elements of the mst
    moment_1_array = np.zeros(Nstar-1)
    moment_2_array = np.zeros(Nstar-1)
    moment_3_array = np.zeros(Nstar-1)
    
    row, column = mst.nonzero()							#the ids of the mst elements which are connected
    for j in range (0, Nstar-1, 1):
	L[j] = mst[row[j], column[j]]
    
    l_mean[i] = np.mean(L)								# average length of the mst 
    
    for j in range (0, Nstar-1, 1):
	moment_1_array[j] = (L[j] - l_mean[i])
	moment_2_array[j] = (L[j] - l_mean[i])**2
	moment_3_array[j] = (L[j] - l_mean[i])**3
	
   
    moment_1[i] = np.mean(moment_1_array)
    moment_2[i] = np.mean(moment_2_array)
    moment_3[i] = np.mean(moment_3_array)
    
  
    
    sigma[i] = math.sqrt(moment_2[i])							
    gamma[i] = moment_3[i]/(sigma[i])**3
    print scipy.stats.skewtest(L)
     
    i = i+1

print "spearmanr test for mean", scipy.stats.spearmanr(l_mean[1:], np.arange(100, 990, 10))
print "spearmanr test for standard deviation", scipy.stats.spearmanr(sigma[1:], np.arange(100, 990, 10)) 
print "spearmanr test for skewness", scipy.stats.spearmanr(gamma[1:], np.arange(100, 990, 10))
'''
plt.plot(np.arange(90, 990, 10), l_mean, markersize=2, color='blue')
plt.ylabel('length mean[pc]', fontsize=15)
plt.xlabel('time[Myr]', fontsize=15)
plt.savefig("../plots/noPE_PI_SN/mean_mst_less40.png")
plt.clf()
plt.plot(np.arange(90, 990, 10), sigma, markersize=2, color='blue')
plt.ylabel('sigma[pc]', fontsize=15)
plt.xlabel('time[Myr]', fontsize=15)
plt.savefig("../plots/noPE_PI_SN/sigma_mst_less40.png")
plt.clf()
plt.plot(np.arange(90, 990, 10), gamma, markersize=2, color='blue')
plt.ylabel('skewness', fontsize=15)
plt.xlabel('time[Myr]', fontsize=15)
plt.savefig("../plots/noPE_PI_SN/skewness_mst_less40.png")
'''

   
      

    
