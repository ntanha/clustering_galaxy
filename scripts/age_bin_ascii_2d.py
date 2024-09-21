# this is a script to read into gadget snapshots and bin the positions of the stars based on their age and write them to an ascii file (2D)

import numpy as np
import readsnap as rs
import math

unit_time_in_Myr = 1e3

num = str(990)

filename_path = "/home/tanha/cheops_mnt/noPE_PI_SN"
 
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

print "number of alive star particles" + str(np.shape(index_alive))
print "length of the x array is" + str(x_star.size)
print "length of the y array is" + str(y_star.size)

x_max = np.amax(x_star)
x_min = np.amin(x_star)
y_max = np.amax(y_star)
y_min = np.amin(y_star)

x_max = str(x_max)
x_min = str(x_min)
y_max = str(y_max)
y_min = str(y_min)

print "boundaries of the snapshot %s are %s %s %s %s"  %(num, x_max, x_min, y_max, y_min)

a = range(0,10,1)
b = range(10,100,10)
c = range(100,1000,100)
d = range(1000,1001,1)

age_star = age_star[index_alive]

tot = 0

# This for loop block goes throuh logarithmic age bins starting from 10 Myrs up to 1000 Myrs

for i in range (0, 18, 1):
  print i
  print "the star particles between age %d and %d Myrs" %((b+c+d)[i], (b+c+d)[i+1])
#  print "the star particles between age %d and %d Myrs" %(0, (a+b+c+d)[i+1])

  id_bin_i = np.where(((b+c+d)[i] <= age_star) & (age_star <= (b+c+d)[i+1]))
 # id_bin_i = np.where((0 <= age_star) & (age_star <= (a+b+c+d)[i+1]))
  age_star_i = age_star[id_bin_i]
  
  print id_bin_i
  print age_star_i

  x_star_i = x_star[id_bin_i]
  y_star_i = y_star[id_bin_i]
  print "The length of x array in the %d bin is %d" %(i, x_star_i.size)
  
  Nstar_i = age_star_i.size
  print "The number of star particles in the %dth bin is %d" %(i, Nstar_i)
  tot = tot + Nstar_i

  Nstar_i = str(Nstar_i)

  arr = np.array((x_star_i,y_star_i),dtype=float)
  
  f = open('./positions/noPE_PI_SN/position_%s_%d_%dto%d_2D.txt' %(num, i,(b+c+d)[i], (b+c+d)[i+1]), 'w')  
  f.write(Nstar_i + '\n')
  f.write(x_max + '\n')
  f.write(x_min + '\n')
  f.write(y_max + '\n')
  f.write(y_min + '\n')
  f.close()

  f = open('./positions/noPE_PI_SN/position_%s_%d_%dto%d_2D.txt' %(num, i,(b+c+d)[i], (b+c+d)[i+1]), 'a') 
  np.savetxt(f,arr)


# this bit does the 0 to 5 and 5 to 10 Myr binning

print "the star particles between age 0 and 10 Myrs"
id_bin_less10 = np.where(age_star <= 10)
age_star_less10 = age_star[id_bin_less10]

x_star_less10 = x_star[id_bin_less10]
y_star_less10 = y_star[id_bin_less10]
print "The length of x array in the first bin is %d" %x_star_less10.size

Nstar_less10 = age_star_less10.size
print "The number of star particles in the first bin is %d" %Nstar_less10
tot = tot + Nstar_less10

Nstar_less10 = str(Nstar_less10)

arr = np.array((x_star_less10, y_star_less10), dtype=float)

f = open('./positions/noPE_PI_SN/position_990_less10.txt', 'w')
f.write(Nstar_less10 + '\n')
f.write(x_max + '\n')
f.write(x_min + '\n')
f.write(y_max + '\n')
f.write(y_min + '\n')
f.close()

f = open('./positions/noPE_PI_SN/position_990_less10.txt', 'a')
np.savetxt(f,arr)

print tot ==  x_star.size

print "the star particles between age 0 and 5 Myrs"
id_bin_less5 = np.where(age_star <= 5)
age_star_less5 = age_star[id_bin_less5]

x_star_less5 = x_star[id_bin_less5]
y_star_less5 = y_star[id_bin_less5]
print "The length of x array in the first bin is %d" %x_star_less5.size

Nstar_less5 = age_star_less5.size
print "The number of star particles in the first bin is %d" %Nstar_less5


Nstar_less5 = str(Nstar_less5)

arr = np.array((x_star_less5, y_star_less5), dtype=float)

f = open('./positions/noPE_PI_SN/position_990_less5.txt', 'w')
f.write(Nstar_less5 + '\n')
f.write(x_max + '\n')
f.write(x_min + '\n')
f.write(y_max + '\n')
f.write(y_min + '\n')
f.close()

f = open('./positions/noPE_PI_SN/position_990_less5.txt', 'a')
np.savetxt(f,arr)


# this part should do the age split binning of 0 to 40 Myrs and bigger than 40 Myrs
tot = 0
print "the star particles between age 40 and 990 Myrs"
id_bin_more40 = np.where(age_star >= 40)
age_star_more40 = age_star[id_bin_more40]

x_star_more40 = x_star[id_bin_more40]
y_star_more40 = y_star[id_bin_more40]
print "The length of x array in the last bin is %d" %x_star_more40.size

Nstar_more40 = age_star_more40.size
print "The number of star particles in the last bin is %d" %Nstar_more40
tot = tot + Nstar_more40

Nstar_more40 = str(Nstar_more40)

arr = np.array((x_star_more40, y_star_more40), dtype=float)

f = open('./positions/noPE_PI_SN/position_990_more40.txt', 'w')
f.write(Nstar_more40 + '\n')
f.write(x_max + '\n')
f.write(x_min + '\n')
f.write(y_max + '\n')
f.write(y_min + '\n')
f.close()

f = open('./positions/noPE_PI_SN/position_990_more40.txt', 'a')
np.savetxt(f,arr)

print "the star particles between age 0 and 40 Myrs"
id_bin_less40 = np.where(age_star <= 40)
age_star_less40 = age_star[id_bin_less40]

x_star_less40 = x_star[id_bin_less40]
y_star_less40 = y_star[id_bin_less40]
print "The length of x array in the last bin is %d" %x_star_less40.size

Nstar_less40 = age_star_less40.size
print "The number of star particles in the last bin is %d" %Nstar_less40
tot = tot + Nstar_less40

Nstar_less40 = str(Nstar_less40)

arr = np.array((x_star_less40, y_star_less40), dtype=float)

f = open('./positions/noPE_PI_SN/position_990_less40.txt', 'w')
f.write(Nstar_less40 + '\n')
f.write(x_max + '\n')
f.write(x_min + '\n')
f.write(y_max + '\n')
f.write(y_min + '\n')
f.close()

f = open('./positions/noPE_PI_SN/position_990_less40.txt', 'a')
np.savetxt(f,arr)
print tot ==  x_star.size
