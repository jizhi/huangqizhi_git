from ArrayAxis import *
from PoolFor import *
from FuncFit import *



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
	n1, n2 = iterable[0]
	array = iterable[1].T
	weight, lw, dnwa, dnwa1, dnwa2, appendsize = iterable[2]
	shapew = np.ones(len(array.shape), int)
	shapew[0] = lw
	weight = weight.reshape(shapew)
	b = array*0.
	#--------------------------------------------------
	def func(x, p) : return p[0]*x + p[1]
	dal = array[:int(len(array)*appendsize)]
	dar = array[-len(dal):]
	pl, pr = [], []
	for j in xrange(array.shape[1]) : 
		pl.append(Leastsq(func, np.arange(len(dal)), dal[:,j], [1,0]))
		pr.append(Leastsq(func, np.arange(len(array)-len(dal), len(array)), dar[:,j], [1,0]))
	dal = np.zeros((lw/2+1,)+array.shape[1:], array.dtype)
	dar = dal*0
	for j in xrange(array.shape[1]) : 
		dal[:,j] = func(np.arange(-len(dal), 0), pl[j])
		dar[:,j] = func(np.arange(len(array), len(array)+len(dal)), pr[j])
	#--------------------------------------------------
	for i in xrange(len(array)) : 
		if (i < lw/2) : 
			dn = lw/2-i
		#	da = np.concatenate([array[:1] for j in xrange(dn)])
			print '2'
			da = dal[-dn:]
			ai = np.concatenate([da, array[:lw/2+1+i]])
		elif (i >= len(array)-lw/2) : 
			dn = lw/2 - (len(array)-1 - i)
		#	da = np.concatenate([array[-1:] for j in xrange(dn)])
			da = dar[:dn]
			ai = np.concatenate([array[i-lw/2:], da])
		else : ai = array[i-lw/2:i+lw/2+1]
		b[i] = (ai * weight).sum(0)
	if (dnwa) : b = b[dnwa1:-dnwa2]
	return b



def Smooth( array, axis, per, times=1, appendsize=0.1, sigma=False, reduceshape=False, Nprocess=1 ) : 
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
	# ND to 2D
	if (len(shape0) == 1) : array = array[:,None]
	elif (len(shape0) > 2) : array = array.reshape(shape0[0], np.prod(shape0[1:]))
	shape1 = array.shape  # 2D, along axis=0
	#--------------------------------------------------
	if (not reduceshape) : 
		if (shape1[1] == 1) : Nprocess = 1
		else : 
			Nprocess = NprocessCPU(Nprocess)[0]
			if (Nprocess > shape1[1]) : Nprocess = shape1[1]
		weight = SmoothWeight( per, times )
		lw = len(weight)
		#--------------------------------------------------
		dnwa, dnwa1, dnwa2 = False, None, None
		if (shape1[0] < lw) : 
			dnwa = True
			dnwa1 = (lw-len(array))/2
			dnwa2 = lw-len(array) - dnwa1
			if (dnwa1 > 0) : a1 = np.concatenate([array[:1]  for i in xrange(dnwa1)])
			else : a1 = array[:0]
			a2 = np.concatenate([array[-1:] for i in xrange(dnwa2)])
			array = np.concatenate([a1, array, a2])
			a1 = a2 = 0 #@
		#--------------------------------------------------
		if (Nprocess == 1) : 
			iterable = ((0,0), array.T, [weight, lw, dnwa, dnwa1, dnwa2, abs(appendsize)])
			b = _Multiprocess_Smooth(iterable)
		else : 
			pool = PoolFor(0, shape1[1], Nprocess)
			b = pool.map_async(_Multiprocess_Smooth, array.T, [weight, lw, dnwa, dnwa1, dnwa2, abs(appendsize)])
			b = np.concatenate(b, 1)
		if (len(shape0) != 2) : b = b.reshape(shape0)
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


