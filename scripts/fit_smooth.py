# run with python fit_smooth.py "sim_type" "age_bin"
# sim_type is one of the following options: PE_PI_SN, noPE_PI_SN, PE_noPI_SN, PE_PI_SN3Myr
# age_bin is one of the options 10Myrs, 100Myrs, 5Myrs, 40Myrs
# this script reads in the ascii file procudec by the 2point2Dlogbins.cpp routine
# The input file is and array containing bins, dd, dr, rr and two point function array in that order
# here we read int the input file, smooth the data points, to make guesses for the fitting. then we use the unsmoothed data to fit (Seamus), then plot.
# we add poisson errors(sqrt(N)) to the data points. Since we are fitting and plotting in logaritmic scale, we need to convert the errors. Based on https://faculty.washington.edu/stuve/log_error.pdf we calculate the arror as d(log(y)) = 0.434(dy/y)

import numpy as np
import matplotlib 
import matplotlib.pyplot as plt
from io import StringIO
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from scipy.optimize import curve_fit
from astropy.convolution import convolve, Gaussian1DKernel
import sys



start = 1                                                       # the index where we start fitting. We choose it to be one (not zero) because all thae data in the first bin are less than resolution in the simulation and not trust worthy
end = 11
tol = 0.5                                                       # tolerance -> controls with how much slope difference the code adds a new power law
smooth = 2                                                      # standard deviation of the smoothing Gaussian kernel

sim_type = str(sys.argv[1])
age_bin = str(sys.argv[2])

def g(x, a, c):
    return (a*x + c)

def f(x, x0, y0, k1, k2):
    return np.piecewise(x, [x<x0, x>x0], [lambda x: k1*(x-x0) + y0, lambda x: k2*(x-x0) + y0])

color = ["dimgray", "indianred", "tomato", "chocolate", "burlywood", "y", "forestgreen", "darkseagreen", "lightskyblue", "cornflowerblue", "darkmagenta", "mediumorchid", "palevioletred", "lightpink", "crimson", "gold", "yellow", "navy", "lightcoral", "firebrick", "darkorange", "tan", "yellowgreen", "darkgreen", "darkcyan", "steelblue", "royalblue", "orchid"]

#Label = [" < 40 Myr", " > 40 Myr"]
Label = [" < 5 Myr", " < 10 Myr"]

filename_base = "/home/nassim/analysis/clustering/2point_data/"+ sim_type +"/2point_990_log_2D_"

l = 0

### defining an array which saves the fitting parameters 
fit_broken = np.zeros((9, 4))
fit_single = np.zeros((9, 2))
error_broken = np.zeros((9, 4))
error_single = np.zeros((9, 4))

### defining an array which keeps track of the reduced xi squared result. The elements are zero if broken power law is a better fit and 1 if the single power law is the better fit
check = np.zeros(9)

plt.clf()
#for k in range(10, 100, 10):
for k in range(0, 2, 1):
    num = str(k) + "to" + str(k+10)
    if k == 0:
#        num = 'less40'
        num = 'less5'
    elif k == 1:
#        num = 'more40'
        num = 'less10'
    filename = filename_base + num + '.txt'
  
    ### read in the input files 
    print ''
    print 'Read data in age range ', num
    
    TwoPoint= np.zeros(12)
    bins = np.arange(12)
    
    data = np.genfromtxt(filename, dtype=float)
    bins = data[:, 0]
    TwoPoint = 1 + data[:, 4]
    
    bins = bins[TwoPoint>10**(-5)]
    TwoPoint = TwoPoint[TwoPoint>10**(-5)]
        
    logbin = np.log10(bins)
    logTwoPoint = np.log10(TwoPoint)
    
    end = np.size(logTwoPoint) -1
    
    ### calculating the error for each data point 
    error = 0.5*(np.log10(TwoPoint + np.sqrt(TwoPoint)) - np.log10(TwoPoint - np.sqrt(TwoPoint)))
    
    ### smooth the data 
    gk = Gaussian1DKernel(smooth)
    logTwoPoint_conv = convolve(logTwoPoint,gk)
    
    ### calculating the derivatives 
    n = np.size(logTwoPoint)
    dy = np.zeros(n-4)                                      # stores the differential between two neigboring points - > n-1 because it's the difference, and we remove the first and the two last points from the data  
    dx = (logbin[n-3] - logbin[1])/n
    g1 = np.zeros(n-5)
    
    for i in range (0, n-4, 1):
        dy[i] = (logTwoPoint_conv[i+1] - logTwoPoint_conv[i])/(dx)
    print 'slopes of the two point correlation function vs. postition bins in log log space are ', dy

    for i in range (0, n-5, 1):
        g1[i] = dy[i]/dy[i+1]
    print 'ratios of the neiboring slopes are ', g1

    
    ### guess the breaking point based on the maximum slope of the line connecting the neigboring points
    g1max = np.max(np.absolute(g1))                                      
    id_guess = np.where(np.absolute(g1) == g1max)                         
    slope1_guess = dy[id_guess[0]+1]
    slope2_guess = dy[id_guess[0]+2]
    bp_x_guess = logbin[id_guess[0]+2]
    bp_y_guess = logTwoPoint[id_guess[0]+2]
    
    slope_guess = (logTwoPoint[end]-logTwoPoint[start])/(logbin[start]-logbin[end])
    c_guess = bp_y_guess - bp_x_guess*slope_guess
    print 'guessed values for the first and second slopes, coordinates of the breaking point and the slope of the linear fits are', slope1_guess, slope2_guess, bp_x_guess, bp_y_guess, slope_guess

    
    ### fitting the broken pawer law
#    popt, pcov = curve_fit(f, logbin[start:end], logTwoPoint[start:end], p0=[bp_x_guess, bp_y_guess, slope1_guess, slope2_guess], sigma = error[start:end])
    popt, pcov = curve_fit(f, logbin[start:end], logTwoPoint[start:end], p0=[bp_x_guess, bp_y_guess, slope1_guess, slope2_guess])
    print '[x0, y0, k1, k2] = ', popt
    fit_broken[l, :] = popt
    for j in range (0, 4, 1):
        error_broken[l, j] = pcov[j, j]
        print 'errors in the ', j, 'th parameter for the broken power law fit is ', pcov[j, j]

        
    ### calculating the reduced chi squared for the broken pawer law
    r_ch_sq_br = np.sum((logTwoPoint[start:end] - f(logbin[start:end], *popt))**2)/((end-start+1)-5)
    
    ### ploting the broken pawer law
    plt.figure(1)
    plt.plot(logbin[start:end], f(logbin[start:end], *popt), '--', markersize=2, c=color[l])
#    plt.plot(logbin[start:end], logTwoPoint[start:end], markersize=4, c=color[l], label="%d Myr < age < %d Myr" %(k, k+10))
#    plt.errorbar(logbin[start:end], logTwoPoint[start:end], yerr = error_broken[l, :], markersize=4, c=color[l], label="%d Myr < age < %d Myr" %(k, k+100))
    plt.plot(logbin[start:end], logTwoPoint[start:end], markersize=4, c=color[l], label=Label[l])

    
    ### fitting the single pawer law
#    popt, pcov = curve_fit(g, logbin[start:end], logTwoPoint[start:end], p0 = [slope_guess, c_guess], sigma = error[start:end])
    popt, pcov = curve_fit(g, logbin[start:end], logTwoPoint[start:end], p0 = [slope_guess, c_guess])
    print '[a, c] = ', popt
    fit_single[l, :] = popt
    for j in range (0, 2, 1):
        error_single[l, j] = pcov[j, j]
        print 'errors in the ', j, 'th parameter for the single pawer law fit is ', pcov[j, j]    
 
 
    ### calculating the reduced chi squared for the single power law
    r_ch_sq_si = np.sum((logTwoPoint[start:end] - g(logbin[start:end], *popt))**2)/((end-start+1)-3)
    
    ### ploting the single pawer law
    plt.figure(2)
    plt.plot(logbin[start:end], g(logbin[start:end], *popt), '--', markersize=2, c=color[l])
#    plt.plot(logbin[start:end], logTwoPoint[start:end], markersize=4, c=color[l], label="%d Myr < age < %d Myr" %(k, k+10))
#    plt.errorbar(logbin[start:end], logTwoPoint[start:end], yerr = error_single[l, :], markersize=4, c=color[l], label="%d Myr < age < %d Myr" %(k, k+100))
    plt.plot(logbin[start:end], logTwoPoint[start:end], markersize=4, c=color[l], label=Label[l])

    
    ### comparing the two fits
    print 'reduced chi squared for the broken power law is ', r_ch_sq_br, 'reduced chi squared for the single power law is ', r_ch_sq_si
    if r_ch_sq_br < r_ch_sq_si :
        check[l] = 0
        print 'broken power law is a better fit'
    else:
        check[l] = 1
        print 'single power law is a better fit'
        
    l = l+1


np.save("../fits/" + sim_type + "/fit_broken_"+ age_bin +"_"+ sim_type +".npy", fit_broken)
np.save("../fits/" + sim_type + "/fit_single_"+ age_bin +"_"+ sim_type +".npy", fit_single)
np.save("../fits/" + sim_type + "/error_broken"+ age_bin +"_"+ sim_type +".npy", error_broken)
np.save("../fits/" + sim_type + "/error_single"+ age_bin +"_"+ sim_type +".npy", error_single)
np.save("../fits/" + sim_type + "/check_"+ age_bin +"_"+ sim_type +".npy", check)

   
plt.figure(1)
plt.legend(loc="best")
plt.ylabel('log Two Point Function + 1', fontsize=15)
plt.xlabel('log distance [pc]', fontsize=15)
plt.savefig("../plots/"+ sim_type +"/twopoint_log_"+ age_bin +"_broken.png")  


plt.figure(2)   
plt.legend(loc="best")
plt.ylabel('log Two Point Function + 1', fontsize=15)
plt.xlabel('log distance [pc]', fontsize=15)
plt.savefig("../plots/"+ sim_type +"/twopoint_log_"+ age_bin +"_single.png")
    
plt.show()

              
      

  
