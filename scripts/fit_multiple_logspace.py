# This script reads in the ascii file produced by 2point2D.cpp routine. 
#The input file is an array containing bins, dd, dr, rr and two point fucntion array in that order. 
#Here we read in the input file, fit a broken linear function to them, and  plot them.

import numpy as np
import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
from io import StringIO
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from scipy.optimize import curve_fit

filename_base = "/home/nassim/analysis/clustering/2point_data/noPE_PI_SN/2point_990_log_2D_"

def f(x, x0, y0, k1, k2):
  return np.piecewise(x, [x < x0], [lambda x:k1*x + y0-k1*x0, lambda x:k2*x + y0-k2*x0])

color = ["dimgray", "indianred", "tomato", "chocolate", "burlywood", "y", "forestgreen", "darkseagreen", "lightskyblue", "cornflowerblue", "darkmagenta", "mediumorchid", "palevioletred", "lightpink", "crimson", "gold", "yellow", "navy", "lightcoral", "firebrick", "darkorange", "tan", "yellowgreen", "darkgreen", "darkcyan", "steelblue", "royalblue", "orchid"]
Label = [" < 40 Myr", " > 40 Myr"]
i = 0
#for k in range(0, 2, 1):
for k in range(10, 20, 10):
#  num = str(k)
  num = str(k) + "to" + str(k+10)
#  if k == 0:
#    num = "less40"
#  elif k == 1:
#    num = "more40"
  
  filename = filename_base + num + '.txt'
  
  print ''
  print 'Read data in age range ', num
  
  TwoPoint= np.zeros(12)
  bins = np.arange(12)
  
  data = np.genfromtxt(filename, dtype=float)
  bins = data[:, 0]
  TwoPoint = 1 + data[:, 4]
  
  bins = bins[TwoPoint>10**(-3)]
  TwoPoint = TwoPoint[TwoPoint>10**(-3)]
    
  logbin = np.log10(bins)
  logTwoPoint = np.log10(TwoPoint)
  
  print logbin
  print logTwoPoint
  
  def f(x, x0, y0, k1, k2):
    return np.piecewise(x, [x < x0], [lambda x:k1*x + y0-k1*x0, lambda x:k2*x + y0-k2*x0])
  
  popt, pcov = curve_fit(f, logbin[1:10], logTwoPoint[1:10])
  print '[x0, y0, k1, k2] = ', popt
  for j in range (0, 4, 1):
    print 'errors in the ', j, 'th parameter is = ', pcov[j, j]
  
  plt.plot(logbin[1:12], logTwoPoint[1:12], markersize=4, c=color[i], label="%d Myr < age < %d Myr" %(k, k+10))
#  plt.plot(logbin[1:12], logTwoPoint[1:12], markersize=4, c=color[i], label=" < %d Myr" %k)
#  plt.plot(logbin[1:12], logTwoPoint[1:12], markersize=4, c=color[i], label=Label[i])
  plt.plot(logbin[1:12], f(logbin[1:12], *popt), '--', markersize=2, c=color[i])

  i = i+1
  
plt.legend(loc="best")
plt.ylabel('log Two Point Function + 1', fontsize=15)
plt.xlabel('log distance [pc]', fontsize=15)
#plt.show()
plt.savefig('../plots/noPE_PI_SN/twopoint_log_10Myrs_piecewisefit.png')  
plt.clf()

i = 0
#for k in range(0, 2, 1):
for k in range(10, 20, 10):
#  num = str(k)
  num = str(k) + "to" + str(k+10)
#  if k == 0:
#    num = "less40"
#  elif k == 1:
#    num = "more40"
  
  filename = filename_base + num + '.txt'
  
  print ''
  print 'Read data in age range ', num
  
  TwoPoint= np.zeros(12)
  bins = np.arange(12)
  
  data = np.genfromtxt(filename, dtype=float)
  bins = data[:, 0]
  TwoPoint = 1 + data[:, 4]

  bins = bins[TwoPoint>10**(-3)]  
  TwoPoint = TwoPoint[TwoPoint>10**(-3)]
  
  logbin = np.log10(bins)
  logTwoPoint = np.log10(TwoPoint)
  
  def g(x, a, b):
      return a*x + b
  popt, pcov = curve_fit(g, logbin[1:10], logTwoPoint[1:10])
  print '[a, b] = ', popt
  for j in range (0, 2, 1):
    print 'errors in the ', j, 'th parameter is = ', pcov[j, j]  
  
  plt.plot(logbin[1:12], logTwoPoint[1:12], markersize=4, c=color[i], label="%d Myr < age < %d Myr" %(k, k+10))
  plt.plot(logbin[:], g(logbin[:], *popt), '--', markersize=2, c=color[i])


  i = i+1
  
plt.legend(loc="best")
plt.ylabel('log Two Point Function + 1', fontsize=15)
plt.xlabel('log distance [pc]', fontsize=15)
#plt.show()
plt.savefig('../plots/noPE_PI_SN/twopoint_log_10Myrs_linearfit_curve.png')  

