# this is a script to read into multiple snapshots and make the SFR vs. time and over plot them on the mst everag/ edge length


import numpy as np
import readsnap as rs
import scipy
from scipy.optimize import curve_fit
from scipy.sparse.csgraph import minimum_spanning_tree
from scipy.stats import skewtest
from scipy.stats import spearmanr
from astropy.convolution import convolve, Gaussian1DKernel 
import math
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm

unit_time_in_Myr = 1e3
unit_mass_in_solarmass = 1e10
unit_time_in_yr = 1e9


begin = 90									#in Myr time of the first snapshot
end = 990									#in Myr time of the last snapshot

sim_type = "PE_PI_SN_3Myr"
    
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
	
    filename_path = "/home/nassim/analysis/clustering/dd_less40/pe_pi_sn_3myr"
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

np.save('/home/nassim/analysis/clustering/mst_stat/' + sim_type + '/l_mean.npy', l_mean)
np.save('/home/nassim/analysis/clustering/mst_stat/' + sim_type + '/sigma.npy', sigma)
np.save('/home/nassim/analysis/clustering/mst_stat/' + sim_type + '/gamma.npy', gamma)

plt.subplot(411)
plt.plot(np.arange(90, 990, 10), l_mean, markersize=1.5, color='blue')
plt.ylabel('length mean[pc]', fontsize=12)
plt.xlabel('time[Myr]', fontsize=12)

plt.subplot(412)
plt.plot(np.arange(90, 990, 10), sigma, markersize=1.5, color='blue')
plt.yticks(np.arange(15, 65, 10))
plt.ylabel('sigma[pc]', fontsize=12)
plt.xlabel('time[Myr]', fontsize=12)

plt.subplot(413)
plt.plot(np.arange(90, 990, 10), gamma, markersize=1.5, color='blue')
plt.yticks(np.arange(1.5, 5.5, 1))
plt.ylabel('skewness', fontsize=12)
plt.xlabel('time[Myr]', fontsize=12)

filename_path = "/home/nassim/cheops_mnt/" + sim_type

filename_base = filename_path + "/snap_MW_mres_nort_"

sfr = np.zeros((end - begin)/10)

i = 0
jump = 1

for k in range (begin, end, 10):
    kk = k*jump
    if (kk < 10):
	num = '00' + str(kk)
    elif (kk < 100):
	num = '0' + str(kk)
    elif (kk < 1000):
	num = str(kk)

    filename = filename_base + num

    print ''
    print 'Read snaphot ', num
    head= rs.snapshot_header(filename)
    
    Ngas   = head.npart[0]    
    Nstar  = head.npart[4]
    
    print "number of gas particle %d" % Ngas
    print "number of star particles %d" % Nstar

    time = np.float64(head.time * unit_time_in_yr)				#in Myr

    print "this snapshot is taken after %f yr" % time   
    
    if (Ngas > 0):
      sfr_arr = rs.read_block(filename, "SFR ", parttype=0) * ((unit_mass_in_solarmass)/(unit_time_in_yr))
      sfr[i] = np.sum(sfr_arr)
      
      i = i+1
      
print "sfr = ", sfr

np.save('/home/nassim/analysis/clustering/mst_stat/' + sim_type + '/sfr.npy', sfr)

plt.subplot(414)
plt.plot(np.arange(begin, end, 10), sfr, markersize=1.5, color='blue')
plt.yticks(np.arange(0, 5e-3, 1e-3))
plt.ylabel('SFR[solar mass/yr]')
plt.xlabel('time[Myr]')

#plt.savefig("../plots/noPE_PI_SN/sfr_overplot.png")
plt.show()
plt.clf()


# here we convolve all the data with a Gaussian1DKernel to smooth the graphs. 
gk = Gaussian1DKernel(2)
l_conv = convolve(l_mean,gk)
sfr_conv = convolve(sfr,gk)
sigma_conv = convolve(sigma,gk)
gamma_conv = convolve(gamma,gk)

# This bit is finding the minima of the convolved mean length of the minimum spanning tree so we can make vertical lines to compare the behaviour of the other paramteres at those moments
t = np.arange(90, 990, 10)
maxima_lconv = np.r_[True, l_conv[1:] < l_conv[:-1]] & np.r_[l_conv[:-1] < l_conv[1:], True]
id_maxima_lconv = np.where(maxima_lconv == True)
print "maximas id = ", id_maxima_lconv
vline = t[id_maxima_lconv]

np.save('/home/nassim/analysis/clustering/mst_stat/' + sim_type + '/l_conv.npy', l_conv)
np.save('/home/nassim/analysis/clustering/mst_stat/' + sim_type + '/sigma_conv.npy', sigma_conv)
np.save('/home/nassim/analysis/clustering/mst_stat/' + sim_type + '/gamma_conv.npy', gamma_conv)
np.save('/home/nassim/analysis/clustering/mst_stat/' + sim_type + '/sfr_conv.npy', sfr_conv)


plt.subplot(411)
plt.plot(np.arange(90, 990, 10), l_conv, markersize=1.5, color='blue')
for l in vline:
  plt.axvline(x=l, color='black')
plt.ylabel('length mean[pc]', fontsize=12)
plt.xlabel('time[Myr]', fontsize=12)

plt.subplot(412)
plt.plot(np.arange(90, 990, 10), sigma_conv, markersize=1.5, color='blue')
plt.yticks(np.arange(15, 65, 10))
for l in vline:
  plt.axvline(x=l, color='black')
plt.ylabel('sigma[pc]', fontsize=12)
plt.xlabel('time[Myr]', fontsize=12)

plt.subplot(413)
plt.plot(np.arange(90, 990, 10), gamma_conv, markersize=1.5, color='blue')
plt.yticks(np.arange(1.5, 5.5, 1))
for l in vline:
  plt.axvline(x=l, color='black')
plt.ylabel('skewness', fontsize=12)
plt.xlabel('time[Myr]', fontsize=12)

plt.subplot(414)
plt.plot(np.arange(begin, end, 10), sfr_conv, markersize=1.5, color='blue')
plt.yticks(np.arange(0, 5e-3, 1e-3))
for l in vline:
  plt.axvline(x=l, color='black')
plt.ylabel('SFR[solar mass/yr]', fontsize=12)
plt.xlabel('time[Myr]', fontsize=12)

#plt.savefig("../plots/noPE_PI_SN/sfr_overplot_convolved.png")    
plt.show()    
    
    
    
    
    
    
    
    
    
    
    
