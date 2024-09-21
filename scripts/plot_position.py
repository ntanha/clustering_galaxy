# this is a script to read into gadget snapshots and bin the positions of the stars based on their age and plot them

import numpy as np
import readsnap as rs
import math
import matplotlib.pyplot as plt

color = ["dimgray", "indianred", "tomato", "chocolate", "burlywood", "y", "forestgreen", "darkseagreen", "lightskyblue", "cornflowerblue", "darkmagenta", "mediumorchid", "palevioletred", "lightpink", "crimson", "gold", "yellow", "navy", "lightcoral", "firebrick", "darkorange", "tan", "yellowgreen", "darkgreen", "darkcyan", "steelblue", "royalblue", "orchid"]

unit_time_in_Myr = 1e3

num = str(990)

filename_path = "/home/nassim/cheops_mnt/noPE_PI_SN"
 
filename_base = filename_path + "/snap_MW_mres_nort_"

filename = filename_base + num

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
age_star = time - SF_time						#in Myr 
print age_star.size

m_imf = rs.read_block(filename, "MIMF", parttype=4)			#a(nstar, 10) shape array. To each star particle onw could assign more than one stellar mass, that's why there is an imf array (1, 10) is assigned to each star particle to show how the mass is distributed between different stars 
index_dead = np.where(m_imf < 0)					#find the entries that are negative (dead stars)
m_imf[index_dead] = 0							#change negative entries with zeros (because there could be still some stellar mass left in imf array of a star particle)
I = np.amax(m_imf, axis=1)						#turn the (nstar, 10) imf arrays into (nstar, 1) arrays that has the maximum imf value
index_alive = np.where(I > 0)						#find the star particles with at least one positive entry in imf array

x_star = x_star[index_alive]
y_star = y_star[index_alive]
age_star = age_star[index_alive]
plt.figure(1, figsize=(10,10))
plt.plot(x_star, y_star, '.', color = 'black')
plt.xlabel('x[pc]')
plt.ylabel('y[pc]')

plt.savefig('/home/nassim/analysis/clustering/plots/compare/star_position_all.png')
plt.show()

tot = 0
"""
k = 0
for i in range (10, 100, 10):
    id_bin_i = np.where((i <= age_star) & (age_star <= i+10))
    age_star_i = age_star[id_bin_i]
    x_star_i = x_star[id_bin_i]
    y_star_i = y_star[id_bin_i]
    Nstar_i = age_star_i.size
    
    tot = tot + Nstar_i
    plt.figure(1, figsize=(10,10))
    plt.plot(x_star_i, y_star_i, '.', color = color[k], label="%d Myr < age < %d Myr" %(i, i+10))
    k = k+1



plt.xlabel('x[pc]')
plt.ylabel('y[pc]')
plt.legend(loc ='best')
plt.savefig('/home/nassim/analysis/clustering/plots/compare/star_position_10.png')
plt.show()




id_bin_10to20 = np.where((10 <= age_star) & (age_star <= 20))
age_star_10to20 = age_star[id_bin_10to20]
x_star_10to20 = x_star[id_bin_10to20]
y_star_10to20 = y_star[id_bin_10to20]
Nstar_10to20 = age_star_10to20.size

id_bin_90to100 = np.where((90 <= age_star) & (age_star <= 100))
age_star_90to100 = age_star[id_bin_90to100]
x_star_90to100 = x_star[id_bin_90to100]
y_star_90to100 = y_star[id_bin_90to100]
Nstar_90to100 = age_star_90to100.size



plt.figure(1, figsize = (13, 6))
plt.clf()
plt.subplot(121)
plt.plot(x_star_10to20, y_star_10to20, '.')
plt.xlabel('x[pc]')
plt.ylabel('y[pc]')
plt.title('10 Myr < age < 20 Myr')
plt.subplot(122)
plt.plot(x_star_90to100, y_star_90to100, '.')
plt.xlabel('x[pc]')
plt.ylabel('y[pc]')
plt.title('90 Myr < age < 100 Myr')

plt.savefig('/home/nassim/analysis/clustering/plots/compare/star_position.png')
"""
#plt.show()

























