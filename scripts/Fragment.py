from pylab import *
import scipy.stats
import numpy
from numpy.fft import rfft
from numpy.fft import fftfreq
from numpy.fft import fftn
from numpy.fft import fftshift
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import minimum_spanning_tree
from scipy.stats import iqr
from scipy.stats import gaussian_kde
from scipy.stats import skew
from scipy.stats import kurtosis
from scipy.stats import ks_2samp
from scipy.stats import anderson_ksamp
#from pyqt_fit import kde
#from pyqt_fit import kde_methods
import ctypes
from numpy.ctypeslib import ndpointer


def AD_test(dist,method,length,num,lower_lim):

	if(method=="NNS"):
		num_cores = len(dist)
	if(method=="MST"):
		num_cores = len(dist) + 1

	numpy.random.seed(1)

	tot_num = 0

	tot_sep = array([],dtype=float)

	while(tot_num<num):

		x = length*numpy.random.random(num_cores)

		y = zeros_like(x)

		pos = column_stack((x,y))

		if(method=="NNS"):
			seps = NNS(pos)
		if(method=="MST"):
			seps, mst = MST(pos)


		if(len(seps)!=sum(seps>lower_lim)):
			continue

		tot_num = tot_num + 1

		tot_sep = concatenate((tot_sep,seps))


	d,c,p = anderson_ksamp([dist,tot_sep])

	return p

def KS_test(dist,method,length,num,lower_lim):

	if(method=="NNS"):
		num_cores = len(dist)
	if(method=="MST"):
		num_cores = len(dist) + 1

	numpy.random.seed(1)

	tot_num = 0

	tot_sep = array([],dtype=float)

	while(tot_num<num):

		x = length*numpy.random.random(num_cores)

		y = zeros_like(x)

		pos = column_stack((x,y))

		if(method=="NNS"):
			seps = NNS(pos)
		if(method=="MST"):
			seps, mst = MST(pos)


		if(len(seps)!=sum(seps>lower_lim)):
			continue

		tot_num = tot_num + 1

		tot_sep = concatenate((tot_sep,seps))

	max_sep = 1.1*amax(tot_sep)
	if(amax(dist) > max_sep):
		max_sep = 1.1*amax(dist)

	'''
	x = linspace(0,max_sep,num_cores)

	kde_dist = gaussian_kde(dist,bw_method="scott")
	k_dist = kde_dist(x)

	kde_tot = gaussian_kde(tot_sep,bw_method="scott")
	k_tot = kde_tot(x)

	c_dist = cumsum(k_dist)/sum(k_dist)
	c_tot = cumsum(k_tot)/sum(k_tot)

	plot(x,c_dist)
	plot(x,c_tot)
	show()
	'''

	a = hist(dist,bins=num_cores,range=(0,max_sep))
	b = hist(tot_sep,bins=num_cores,range=(0,max_sep))

	c_dist = cumsum(a[0])/sum(a[0])
	c_tot = cumsum(b[0])/sum(b[0])

	d,p = ks_2samp(dist,tot_sep)

	return p


def bimodal_coeff(seps):

	g = skew(seps)
	k = kurtosis(seps,fisher=False) - 3

	n = len(seps)

	n_factor = (3*(n-1)**2) / ((n-2)*(n-3))

	beta = (g**2 + 1.0) / (k+n_factor)

	return beta

def Stat_Sig(dist, method, length, num, lower_lim):

	med_sep = []
	mea_sep = []
	std_sep = []
	iqr_sep = []

	beta = []

	if(method=="NNS"):
		num_cores = len(dist)
	if(method=="MST"):
		num_cores = len(dist) + 1

	numpy.random.seed(1)

	tot_num = 0

	while(tot_num<num):


		x = length*numpy.random.random(num_cores)

		y = zeros_like(x)

		pos = column_stack((x,y))

		if(method=="NNS"):
			seps = NNS(pos)
		if(method=="MST"):
			seps, mst = MST(pos)


		if(len(seps)!=sum(seps>lower_lim)):
			continue

		tot_num = tot_num + 1

		med_sep.append(median(seps))
		iqr_sep.append(iqr(seps))

		mea_sep.append(mean(seps))
		std_sep.append(std(seps))

		beta.append(bimodal_coeff(seps))

	med_sep = array(med_sep,dtype=float)
	mea_sep = array(mea_sep,dtype=float)

	iqr_sep = array(iqr_sep,dtype=float)
	std_sep = array(std_sep,dtype=float)

	beta = array(beta,dtype=float)


	med_dist = median(dist)
	iqr_dist = iqr(dist)

	mea_dist = mean(dist)
	std_dist = std(dist)

	beta_dist = bimodal_coeff(dist)



	med_iqr = column_stack((med_sep,iqr_sep))
	mea_std = column_stack((mea_sep,std_sep))


	maxx = amax(med_sep)
	if(med_dist>maxx):
		maxx = med_dist

	maxy = amax(iqr_sep)
	if(iqr_dist>maxy):
		maxy = iqr_dist
	
	kde1 = gaussian_kde1(med_iqr.T,bw_method="scott")
	x1,y1 = mgrid[0:1.1*maxx:100j,0:1.1*maxy:100j]
	po = vstack([x1.ravel(),y1.ravel()])

	k1 = reshape(kde1(po).T,x1.shape)
	k1 = k1.T

	maxx = amax(x1)
	minx = amin(x1)
	dx = (maxx - minx)/ 100.
	
	maxy = amax(y1)
	miny = amin(y1)
	dy = (maxy - miny)/ 100.

	ix = int((med_dist-minx)/dx)
	iy = int((iqr_dist-miny)/dy)

	print minx, maxx, miny, maxy
	print med_dist,iqr_dist
	print ix,iy


	val = k1[iy,ix]
	less = k1[k1<=val]
	'''
	imshow(k1/(numpy.sum(k1)*dx*dy),origin=0,extent=[0,amax(x1),0,amax(y1)])
	colorbar()
	contour(k1/(numpy.sum(k1)*dx*dy),extent=[0,amax(x1),0,amax(y1)],levels=[val/(numpy.sum(k1)*dx*dy)],colors="white")
	plot(med_dist,iqr_dist,"ro",markersize=3)
	show()
	'''
	p_med = numpy.sum(less) / numpy.sum(k1)





	maxx = amax(mea_sep)
	if(mea_dist>maxx):
		maxx = mea_dist

	maxy = amax(std_sep)
	if(std_dist>maxy):
		maxy = std_dist

	kde2 = gaussian_kde(mea_std.T,bw_method="scott")
	x2,y2 = mgrid[0:1.1*maxx:100j,0:1.1*maxy:100j]
	po2 = vstack([x2.ravel(),y2.ravel()])

	k2 = reshape(kde2(po2).T,x2.shape)
	k2 = k2.T

	maxx = amax(x2)
	minx = amin(x2)
	dx = (maxx - minx)/ 100.
	
	maxy = amax(y2)
	miny = amin(y2)
	dy = (maxy - miny)/ 100.

	ix = int((mea_dist-minx)/dx)
	iy = int((std_dist-miny)/dy)

	val = k2[iy,ix]
	less = k2[k2<=val]
	'''
	imshow(k2/(numpy.sum(k2)*dx*dy),origin=0,extent=[0,amax(x2),0,amax(y2)])
	colorbar()
	contour(k2/(numpy.sum(k2)*dx*dy),extent=[0,amax(x2),0,amax(y2)],levels=[val/(numpy.sum(k2)*dx*dy)],colors="white")
	plot(mea_dist,std_dist,"ro",markersize=3)
	show()
	'''
	p_mea = numpy.sum(less) / numpy.sum(k2)


	kde3 = gaussian_kde(beta,bw_method="scott")
	bx = linspace(0,1,10000)
	k3 = kde3(bx)

	dx = bx[1] - bx[0]

	ix = int(beta_dist / dx)

	val = k3[ix]
	less = k3[k3<=val]

	p_beta = numpy.sum(less)/ numpy.sum(k3)


	''''
	plot(bx,k3/(sum(k3)*dx))
	axvline(x=beta_dist,color="k")
	show()
	'''

	

	return p_med, p_mea, p_beta



def NNS(pos):

	n_pos = len(pos[:,0])

	seps = zeros(n_pos)

	for ii in range(0,n_pos):

		min_sep = 1e99
		for jj in range(0,n_pos):

			if(ii==jj):
				continue

			sep = sqrt( (pos[ii,0]-pos[jj,0])**2 + (pos[jj,1]-pos[ii,1])**2 )

			if(sep<min_sep):
				min_sep = sep

		seps[ii] = min_sep

	return seps


def NSS(pos):

	n_pos = len(pos[:,0])

	sep = zeros((n_pos,n_pos))

	for ii in range(0,n_pos):

		for jj in range(0,n_pos):

			sep[ii,jj] = sqrt( (pos[ii,0]-pos[jj,0])**2 + (pos[jj,1]-pos[ii,1])**2 )

		sep[ii,:] = sort(sep[ii,:])

	mean_seps = mean(sep,axis=0)
	mean_std = std(sep,axis=0)

	med = median(sep,axis=0)
	iqr = scipy.stats.iqr(sep,axis=0)

	return sep, mean_seps[1:], mean_std[1:], med[1:], iqr[1:]


def TwoPoint(pos,nruns,bounds):

	n_pos = len(pos[:,0])

	pos_x = pos[:,0]
	pos_y = pos[:,1]

	pos_x = array(pos_x,dtype=double)
	pos_y = array(pos_y,dtype=double)

	x_size=1000

	DD = zeros(x_size, dtype=double)
	DR = zeros(x_size, dtype=double)
	RR = zeros(x_size, dtype=double)

	lib = ctypes.cdll.LoadLibrary("./adder.so")

	twopoint = lib.TwoPoint

	twopoint.restype = None

	twopoint.argtypes = [ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"), ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"), ctypes.c_int, ctypes.c_int, ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"), ctypes.c_int, ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),ndpointer(ctypes.c_double, flags="C_CONTIGUOUS")]

	xmin = bounds[0]
	xmax = bounds[1]
	
	ymin = bounds[2]
	ymax = bounds[3]

	max_dist = sqrt((xmax-xmin)**2+(ymax-ymin)**2)

	twopoint(pos_x,pos_y,n_pos,nruns,bounds,x_size,DD,DR,RR)
	
	w = (DD - 2*DR + RR)/RR	

	sep = linspace(0,max_dist,x_size)

	return sep, w




def ApproxTwoPoint(pos,nruns,bounds,error):

	n_pos = len(pos[:,0])

	pos_x = pos[:,0]
	pos_y = pos[:,1]

	pos_x = array(pos_x,dtype=double)
	pos_y = array(pos_y,dtype=double)

	xmin = bounds[0]
	xmax = bounds[1]
	
	ymin = bounds[2]
	ymax = bounds[3]

	max_dist = sqrt((xmax-xmin)**2+(ymax-ymin)**2)

        xstart = 2.0
	xend = max_dist
	
	x_size = 50
	'''
	Lmin = 2.0

	bins = int((np.log10((xend)/(xstart)))/(np.log10(1+(Lmin)/(xstart)))) + 1
	x_size = bins
	Lmin = xstart*(np.power((xend/xstart),(1.0/(bins-1))) - 1)
	logbinsize = (np.log10(xend) - np.log10(xstart))/(bins-1)
	'''
	
	DD = zeros(x_size, dtype=double)
	DR = zeros(x_size, dtype=double)
	RR = zeros(x_size, dtype=double)

	lib = ctypes.cdll.LoadLibrary("./adder.so")

	twopoint = lib.ApproxTwoPoint

	twopoint.restype = None

	twopoint.argtypes = [ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"), ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"), ctypes.c_int, ctypes.c_int, ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"), ctypes.c_int, ctypes.c_double, ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),ndpointer(ctypes.c_double, flags="C_CONTIGUOUS")]

	twopoint(pos_x,pos_y,n_pos,nruns,bounds,x_size,error,DD,DR,RR)
	
	w = (DD - 2*DR + RR)/RR	
	
        sep = linspace(np.log10(xstart),np.log10(xend),x_size)
        
        '''
	sep = np.zeros(bins)
	
	for i in range(0, bins, 1):
#                sep[i] = np.power(10, (np.log10(xstart) + i*logbinsize))
                sep[i] = np.log10(xstart) + i*logbinsize
	'''
	return sep, w





def MST(pos):

	n_pos = len(pos[:,0])
	
	DD = zeros((n_pos,n_pos))

	for jj in range(0,n_pos):
		for kk in range(jj,n_pos):

			DD[jj,kk] = sqrt( (pos[jj,0]-pos[kk,0])**2 + (pos[jj,1]-pos[kk,1])**2 )
			DD[kk,jj] = DD[jj,kk]



	DD2 = csr_matrix(DD)

	t = minimum_spanning_tree(DD2)

	mst = t.toarray().astype(float)

	r = mst.flatten()
	r = r[r>0]

	dis=[]

	for jj in range(0,len(r)):

		dis.append(r[jj])

	return array(dis,dtype=float), mst

def FT_Map(Map):

	fmap = fftn(Map)

	fmap = fftshift(fmap)

	k = linspace(0,200,201)
	power = zeros_like(k)
	n = zeros_like(k)

	for ii in range(0,len(Map[:,0])):
		for jj in range(0,len(Map[0,:])):

			k2 = sqrt( (ii-200)**2 + (jj-200)**2 )

			k3 = int(k2)

			if(k3<len(k)):

				power[k3] = power[k3] + abs(fmap[ii,jj])
				n[k3] = n[k3] + 1.


	n[n==0] = 1
	power = power/n


	std = zeros_like(k)

	for ii in range(0,len(Map[:,0])):
		for jj in range(0,len(Map[0,:])):

			k2 = sqrt( (ii-200)**2 + (jj-200)**2 )

			k3 = int(k2)

			if(k3<len(k)):

				std[k3] = std[k3] + (abs(fmap[ii,jj]) - power[k3])**2

	std = sqrt(std/(n*(n-1)))

	return k,power,std


def FT_Spine(spine):

	spine = spine - mean(spine)

	fspine = rfft(spine)
	#fspine = fftshift(fspine)

	print len(fspine)

	n = len(fspine)

	k = linspace(0,len(fspine)-1,len(fspine))

	fs = abs(fspine)
	
	return k,fs

	


def DSF(Map, upper_dist, dx, cut):

	cut_index = where(Map>cut)

	num = len(cut_index[0][:])

	xnum = int(upper_dist/dx) + 1

	dsf = zeros((xnum,num))

	for ii in range(0,num-1):

		if(ii%100 == 0):
			print ii, num

		for jj in range(ii+1,num):

			distance = sqrt( (cut_index[0][ii]-cut_index[0][jj])**2 + (cut_index[1][ii]-cut_index[1][jj])**2 )

			d = int(distance)

			if(d<xnum):

				dsf[d,ii] = fabs(Map[cut_index[0][ii],cut_index[1][ii]] - Map[cut_index[0][jj],cut_index[1][jj]])


	DSF = zeros(xnum)
	IQR = zeros(xnum)

	for ii in range(1,xnum):

		t_dsf = array(dsf[ii,:])
		t_dsf = t_dsf[t_dsf>0]

		DSF[ii] = median(t_dsf)

		if(DSF[ii]>0):
			IQR[ii] = float(0.5*( percentile(t_dsf,75) - percentile(t_dsf,25)))

	length = linspace(0,upper_dist, xnum)	

	return length, DSF, IQR

