import matplotlib.pyplot as plt
import numpy as np
import Fragment
from scipy.stats import iqr
from astropy.io import fits
from scipy.optimize import curve_fit
from astrodendro import Dendrogram
from astrodendro.analysis import PPStatistic
import sys
import time

# this bit reads the input out of the folder tree /positions/ and fill the position array pos.
print "reading the input file"
input_file = open('/home/nassim/analysis/clustering/positions/PE_noPI_SN/position_990_16_70to80_2D.txt', 'r')
Nstar = int(input_file.readline())
x_max = float(input_file.readline())
x_min = float(input_file.readline())
y_max = float(input_file.readline())
y_min = float(input_file.readline())
x_array = input_file.readline()
y_array = input_file.readline()
x = x_array.split()
y = y_array.split()
pos = np.zeros((Nstar,2))
for i in range(0, Nstar, 1):
    pos[i, 0] = float(x[i])
    pos[i, 1] = float(y[i])

# ApproxTwoPoint is a function from Fragment.py, which calculates 2 point correlation function with KDE, it calls the functions from add.c So first the c file should be compiled with "gcc -shared -Wl,-soname,adder -o adder.so -fPIC add.c -g". The ApproxTwoPoint gets (pos, nruns, bounds, error) as input and returns (sep, w).  pos is the position array. nruns is the number of random field. bounds is an array with 4 elements(xmin, xmax, ymin, ymax). Error is the opening criteria of the tree (opening distance) and calculated as error > 10000/(Nstar**2 * nruns)

nruns = 1000
bounds = np.zeros(4)
bounds[0] = x_min
bounds[1] = x_max
bounds[2] = y_min
bounds[3] = y_max
error = (10000./((Nstar**2)*nruns))*2

sep, w = Fragment.ApproxTwoPoint(pos, nruns, bounds, error)
print sep
print w

plt.semilogy(sep, w)
plt.show()