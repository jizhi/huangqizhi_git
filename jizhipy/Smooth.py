from ArrayAxis import *
from PoolFor import *



def SmoothWeight( per, times ) : 
	if (per %2 == 0) : 
		print 'Warning: SmoothWeight(), per must be odd, reset per='+str(per)+' to per='+str(per+1)
		per +=1
	a0 = np.ones(per)
	for n in xrange(1, times+1) : 
		if (times == 1) : break
		if (n == 1) : continue
#		a = np.zeros([per, 1+2*n])          #-- Method 1 --
		a = np.zeros([per, 1+2*n*(per/2)])  #-- Method 2 --
		for j in xrange(len(a)) : 
			a[j,j:j+len(a0)] = a0
#		a0 = a.sum(0)
		a0 = a.sum(0) /per
#	weight = a0 / (1.*per)**times
	weight = a0/per  #@ Remember this !!!
	return weight



def _Multiprocess_Smooth( iterable ) : 
	array, weight, lw, dnwa, dnwa1, dnwa2 = iterable
	b = array*0
	for i in xrange(len(array)) : 
		if (i < lw/2) : 
			da = np.concatenate([array[:1] for j in xrange(lw/2-i)])
			ai = np.concatenate([da, array[:lw/2+1+i]])
		elif (i >= len(array)-lw/2) : 
			dn = lw/2 - (len(array)-1 - i)
			da = np.concatenate([array[-1:] for j in xrange(dn)])
			ai = np.concatenate([array[i-lw/2:], da])
		else : ai = array[i-lw/2:i+lw/2+1]
		b[i] = (ai * weight).sum(0)
	if (dnwa) : b = b[dnwa1:-dnwa2]
	return b



def Smooth( array, axis, per, times=1, sigma=False, reduceshape=False, Nprocess=1 ) : 
	'''
	Smooth/Average/Mean array along one axis.
	We can also use spsn.convolve() to do this, but spsn.convolve() will cost much more memory and time, so, the function written here is the best and fastest.

	Weighted average, sigma will reduce to 1/sqrt{per**times}
	Equivalent to average over per**times

	axis:
		array will be smoothed/averaged along which axis.

	per:
		How many bins/elements to average.

	times:
		How many times to smooth.
		For random noise, times=4 is OK
		Note that this parameter is just for reduceshape=False

	sigma:
		False, True, int, np.array
		If False, don't return the error.
		If True, calculate the error from input array.
		If int or np.array, use this sigma to calculate the error of the result.

	reduceshape:
		False, return.shape = array.shape
		True, return.shape < array.shape

	# Also, we can use Convolve to do it:
	#    w = array*0
	#    w = w[len(array)/2-per/2:len(array)/2+per/2+1] = 1
	#    return Convolve(array, w/w.sum())

	Note that:
		Wide of 2 windows: w1 < w2
		a1  = Convolve(a,  w1)
		a2  = Convolve(a,  w2)
		a12 = Convolve(a1, w2)
		=> a12 = Smooth(a2)
		But the beam sizes of the result maps are similar (roughly the same), Beamsize(a12) >= Beamsize(a2).
	'''
	if (per%2 == 0) : per +=1
	if (per<=1 or times<=0) : return array
	array = np.array(array)
	per, atype, shape = int(round(per)), array.dtype, array.shape
	if (axis < 0) : axis = len(shape) + axis
	if (axis >= len(shape)) : Raise(Exception, 'axis='+str(axis)+', array.shape='+str(shape)+', axis out of array.shape')
	# Move axis to axis=0
	array = ArrayAxis(array, axis, 0, 'move')
	shape0, shape = array.shape, array.shape
	# Reduce to 2D
	if (len(shape0) > 2) : 
		array = array.reshape(shape[0], np.prod(shape[1:]))
		shape = array.shape
	#--------------------------------------------------
	if (not reduceshape) : 
		weight = SmoothWeight( per, times )
		lw = len(weight)
		shapew = npfmt(shape)
		shapew[0] = lw
		shapew[1:] = 1
		weight = weight.reshape(shapew)
		#--------------------------------------------------
		dnwa, dnwa1, dnwa2 = False, None, None
		if (len(array) < lw) : 
			dnwa = True
			dnwa1 = (lw-len(array))/2
			dnwa2 = lw-len(array) - dnwa1
			if (dnwa1 > 0) : a1 = np.concatenate([array[:1]  for i in xrange(dnwa1)])
			else : a1 = []
			a2 = np.concatenate([array[-1:] for i in xrange(dnwa2)])
			array = np.concatenate([a1, array, a2])
			a1 = a2 = 0 #@
		#--------------------------------------------------
		Nprocess = NprocessCPU(Nprocess)[0]
		if (Nprocess > shape[1]) : Nprocess = shape[1]
		if (Nprocess == 1) : 
			iterable = (array, weight, lw, dnwa, dnwa1, dnwa2)
			b = _Multiprocess_Smooth(iterable)
		else : 
			nlist = np.linspace(0, shape[1], Nprocess+1).astype(int)
			iterable = []
			for i in xrange(len(nlist)-1) : 
				iterable.append( (array[:,nlist[i]:nlist[i+1]], weight, lw, dnwa, dnwa1, dnwa2) )
			pool = multiprocessing.Pool(Nprocess)
			b = pool.map_async(_Multiprocess_Smooth, iterable).get(10**10)
			b = np.concatenate(b, 1)
		if (len(shape0) > 2) : b = b.reshape(shape0)
	#--------------------------------------------------
	else : 
		if (times == 1) : Npix = len(array)/per
		else : Npix = times
		n = np.linspace(0, len(array), Npix+1).astype(int)
		b = np.zeros((Npix-1,) + shape[1:])
		for i in xrange(len(n)-1) : 
			b[i] = array[n[i]:n[i+1]].mean(0)
	#--------------------------------------------------
	b = ArrayAxis(b, 0, axis, 'move')
	return b


