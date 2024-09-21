# this script reads into multiple dd-array ascii files and produces a minimum spanning tree and a histogram plot for each file and in the end an average (mean and median) histogram of the minimum spanning tree

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from io import StringIO
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from scipy.optimize import curve_fit
from scipy.sparse.csgraph import minimum_spanning_tree
import readsnap as rs

unit_time_in_Myr = 1e3

i = 0
jump = 1
bins = 50
low = np.log10(2.)
high = np.log10(1000.)

beg = 90
end = 990
hist_tot = np.zeros(bins)
hist_array = np.zeros((((end - beg)/10), bins))


for k in range (beg, end, 10):
    kk = k*jump
    if (kk < 10):
	num = '00' + str(kk)
    elif (kk < 100):
	num = '0' + str(kk)
    elif (kk < 1000):
	num = str(kk)
	
    filename_path = "/home/nassim/cheops_mnt/noPE_PI_SN"
 
    filename_base = filename_path + "/snap_MW_mres_nort_"

    filename = filename_base + num
    plt.clf()
    
    print ''
    print 'Read snaphot ', num
    head= rs.snapshot_header(filename)

    Ngas   = head.npart[0]    
    Nstar  = head.npart[4]

    print "number of gas particle %d" % Ngas
    print "number of star particles %d" % Nstar

    time = np.float64(head.time * unit_time_in_Myr)				#in Myr

    print "this snapshot is taken after %f Myr" % time

    if (Nstar > 0):
	pos_star = rs.read_block(filename, "POS ", parttype=4) 
	x_star   = np.float128(pos_star[:,0])*1000.				#in pc
	y_star   = np.float128(pos_star[:,1])*1000.

    SF_time = rs.read_block(filename, "AGE ", parttype=4)			#time of birth
    SF_time = SF_time * unit_time_in_Myr					#in Myr
    age_star = time - SF_time							#in Myr 

    m_imf = rs.read_block(filename, "MIMF", parttype=4)				#a(nstar, 10) shape array. To each star particle onw could assign more than one stellar mass, that's why there is an imf array (1, 10) is assigned to each star particle to show how the mass is distributed between different stars 
    index_dead = np.where(m_imf < 0)						#find the entries that are negative (dead stars)
    m_imf[index_dead] = 0							#change negative entries with zeros (because there could be still some stellar mass left in imf array of a star particle)
    I = np.amax(m_imf, axis=1)							#turn the (nstar, 10) imf arrays into (nstar, 1) arrays that has the maximum imf value
    index_alive = np.where(I > 0)						#find the star particles with at least one positive entry in imf array

    x_star = x_star[index_alive]
    y_star = y_star[index_alive]

    print "number of alive star particles" + str(np.shape(index_alive))
    print "length of the x array is" + str(x_star.size)
    print "length of the y array is" + str(y_star.size)

    age_star = age_star[index_alive]

    id_bin_less40 = np.where (age_star <= 40)
    age_star_less40 = age_star[id_bin_less40]

    x_star_less40 = x_star[id_bin_less40]
    y_star_less40 = y_star[id_bin_less40]
    
    filename_path = "/home/nassim/analysis/clustering/dd_less40/nope_pi_sn"
    filename_base = filename_path + "/dd_less40_"
    filename = filename_base + num + ".txt"

    print ''
    print 'Read dd-array file ', filename
    
    FullTree = np.genfromtxt(filename, dtype=float)
    mst = minimum_spanning_tree(FullTree)
    
    
    Nstar = np.shape(mst)[0]
    print "number of stars in this snapshot between age 0 and 40 Myrs is %d" %Nstar
    
    
    """this bit makes minimum spanning tree plots for each snapshot"""
    print "making the mst plot for this snapshot"
    L = np.zeros(Nstar-1)								#an array of all the non-zero elements of the mst
    row, column = mst.nonzero()							#the id of the elements of the mst which are connected
    for ii in range(0, Nstar-1, 1):
	x_ii = (x_star_less40[row[ii]], x_star_less40[column[ii]])
	y_ii = (y_star_less40[row[ii]], y_star_less40[column[ii]])
	plt.plot(x_ii, y_ii, c = 'blue')
	axes = plt.gca()	
	axes.set_xlim([-2000,2500])
	axes.set_ylim([-2500,2500])
    plt.title('minimum spanning tree for stars of age < 40 myrs')
    plt.savefig('../plots/noPE_PI_SN/mst_%s.png' %num)
    plt.clf()
    
    """this bit makes histogram for each snapshot"""
    print "making histogram for this snapshot"
    for j in range (0, Nstar-1, 1):
	L[j] = mst[row[j], column[j]]
      
    L = np.log10(L[np.where(L!=0)])
    hist_j, bin_j, patch = plt.hist(L, bins,range=(low, high), normed=True)
    plt.xlabel('seperation [log pc]', fontsize=15)
    plt.ylabel('normalized frequency', fontsize=15)
    plt.title('histogram of minimum spanning tree for stars of age < 40 myrs')
    axes = plt.gca()
    axes.set_ylim([0,1])
    plt.savefig('../plots/noPE_PI_SN/histogram_mst_norm_less40_%s.png' %num)
    hist_array[i] =  hist_j
    hist_tot = hist_tot + hist_j
    i = i+1
    
plt.clf()
hist_mean = np.mean(hist_array, axis=0)
print "the average histogram is"
print hist_mean
np.save('../histograms/noPE_PI_SN/average_hist_less40_90to990.npy', hist_mean)
np.save('../histograms/noPE_PI_SN/bins.npy', bin_j)
width = (high - low)/bins
plt.bar(bin_j[0:-1], hist_mean, width,  align='edge')
plt.xlabel('seperation [log pc]', fontsize=15)
plt.ylabel('normalized frequency', fontsize=15)
plt.title('average histogram of minimum spanning tree for stars of age < 40 myrs')
axes = plt.gca()
axes.set_ylim([0,1])
plt.savefig("../plots/noPE_PI_SN/hist_average_less40_90to990.png")
plt.clf()

plt.clf()
hist_median = np.median(hist_array, axis=0)
print "the median histogram is"
print hist_median
np.save('../histograms/noPE_PI_SN/median_hist_less40_90to990.npy', hist_median)
plt.bar(bin_j[0:-1], hist_median, width, align='edge')
plt.xlabel('seperation [log pc]', fontsize=15)
plt.ylabel('normalized frequency', fontsize=15)
plt.title('median histogram of minimum spanning tree for stars of age < 40 myrs')
axes = plt.gca()
axes.set_ylim([0,1])
plt.savefig("../plots/noPE_PI_SN/hist_median_less40_90t0990.png")
    
    
    
    
    
    
    
    
