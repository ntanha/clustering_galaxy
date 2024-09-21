# this is a script to read into multiple gadget snapshots and bin the positions of the stars based on their age and write each of them to ascii files
# the output are the positions of the stars of the same age range throughout the whole simulation, one file per snapshot.

import numpy as np
import readsnap as rs
import math

unit_time_in_Myr = 1e3

filename_path = "/home/nassim/cheops_mnt/noPE_PI_SN"

filename_base = filename_path + "/snap_MW_mres_nort_"

i = 0
jump = 1

for k in range (70, 990, 10):
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

    time = np.float64(head.time * unit_time_in_Myr)				#in Myr

    print "this snapshot is taken after %f Myr" % time

    if (Nstar > 0):
	pos_star = rs.read_block(filename, "POS ", parttype=4) 
	x_star   = np.float128(pos_star[:,0])*1000.				#in pc
	y_star   = np.float128(pos_star[:,1])*1000.

    SF_time = rs.read_block(filename, "AGE ", parttype=4)			#time of birth
    SF_time = SF_time * unit_time_in_Myr					#in Myr
    age_star = time - SF_time							#in Myr 
    print age_star.size

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
    
    print "the star particles between age 0 and 40 Myrs"
    id_bin_less40 = np.where(age_star <= 40)
    age_star_less40 = age_star[id_bin_less40]

    x_star_less40 = x_star[id_bin_less40]
    y_star_less40 = y_star[id_bin_less40]
    print "The length of x array in the first bin is %d" %x_star_less40.size

    Nstar_less40 = age_star_less40.size
    print "The number of star particles in the first bin is %d" %Nstar_less40

    Nstar_less40 = str(Nstar_less40)

    arr = np.array((x_star_less40, y_star_less40), dtype=float)

    f = open('../position_less40/nope_pi_sn/position_less40_%s.txt' %num, 'w')
    f.write(Nstar_less40 + '\n')
    f.write(x_max + '\n')
    f.write(x_min + '\n')
    f.write(y_max + '\n')
    f.write(y_min + '\n')
    f.close()

    f = open('../position_less40/nope_pi_sn/position_less40_%s.txt' %num, 'a')
    np.savetxt(f,arr)
    
    i = i+1
