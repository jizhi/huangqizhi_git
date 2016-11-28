from ArrayAxis import *
from PoolFor import *
#from FuncFit import *
from Gaussian import *
from IsType import *



def _GaussianFilter( per, times ) : 
	'''
	a(n) = S(n) - S(n-1)
	S(n) = a(n) + a(n-1) + a(n-2) + ... + a(2) + a(1)

	per:
		per >=2, can be odd or even

	times:
		times >=1

	return:
		gaussianfilter.size = (per-1) * times + 1
	'''
	per, times = np.array([per, times]).round().astype(int)
	if (per<2 or times<1) : 
		Raise(Warning, 'GaussianFilter(per, times), per='+str(per)+', times='+str(times)+', return gaussianfilter=[1]')
		return np.array([1])
#	a0 = np.ones(per)             # Method 1
	a0 = np.ones(per) / per       # Method 2
	for n in xrange(2, times+1) : 
		n1 = (per-1) * n + 1
		a1 = np.zeros([per, n1])
		for i in xrange(per) : 
			a1[i,i:i+a0.size] = a0
#		a0 = a1.sum(0)             # Method 1
		a0 = a1.sum(0) / per       # Method 2
#	gaussianfilter = a0 / (1.per)**times  # Method 1
	gaussianfilter = a0                   # Method 2
	std = 1. / ((2*np.pi)**0.5 * gaussianfilter.max())
	return [gaussianfilter, std]





def GaussianFilter( per=None, times=None, std=None, Npix=None ) : 
	'''
	(1) per, times
		Use per, times to generate gaussian filter
	(2) sigma, Npix
		Use sigma, npix to generate gaussian filter
	'''
	if (per is not None or times is not None) : 
		gaussianfilter, std = _GaussianFilter(per, times)
	else : 
		x = np.arange(Npix)
		x = x - x.mean()
		gaussianfilter = GaussianValue(x, 0, std)
	return [gaussianfilter, std]





def _Multiprocess_Smooth( iterable ) : 
	array = iterable[1].T  # smooth along axis-0
	weight, nla, nra, sigma = iterable[2]
	weight = weight[:,None]
	n1, n2, n = nla, len(array)-nra, len(weight)
	arr = array*0
	if (sigma is True) : arrstd = array*0
	for i in xrange(n1, n2) : 
		a = array[i-nla:i-nla+n]
		arr[i] = (a * weight).sum(0)
		if (sigma is True): arrstd[i] = ((a-arr[i])*weight).std(0)
	arr = arr[n1:n2]
	if (sigma is True) : 
		arrstd = arrstd[n1:n2]
		arr = np.concatenate([arr, arrstd], 0)
	return arr





def _Multiprocess_SmoothReduce( iterable ) : 
	n, array = iterable[1]
	sigma = iterable[2]
	n = n - n.min()
	arr = np.zeros([len(n), array.shape[1]])
	if (sigma is True) : arrstd = arr*0
	for i in xrange(len(n)) : 
		a = array[n[i,0]:n[i,1]]
		arr[i] = a.mean(0)
		if (sigma is True) : arrstd[i] = a.std(0)
	if (sigma is True) : arr = np.concatenate([arr, arrstd], 1)
	return arr





def Smooth( array, axis, perORfilter, times=1, sigma=False, reduceshape=False, Nprocess=1 ) : 
	'''
	Smooth/Average/Mean array along one axis.
	We can also use spsn.convolve() to do this, but spsn.convolve() will cost much more memory and time, so, the function written here is the best and fastest.

	axis:
		array will be smoothed/averaged along this axis.

	perORfilter:
		(1) per: How many bins/elements to be averaged to get one output point.
		(2) filter: Use this filter directly

	times:
		How many times to smooth.

	sigma:
		False, True, int, np.array
		If False, don't return the error.
		If True, calculate the error of the output.
		If int or np.array, use this sigma to calculate the error of the output.

	reduceshape:
		False, return.shape = array.shape
		True , return.shape < array.shape

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
	if (IsType().isnum(perORfilter)) : 
		which = 'per'
		per = int(round(perORfilter))
		if (per <= 1) : return array
	else : 
		which = 'filter'
		weight = npfmt(perORfilter).flatten()
		if (weight.size == 1) : return array
	times = int(round(times))
	if (times <= 0) : return array
	# per is even, b[i] = a[i+1-per/2:i+1+per/2]
	# per is  odd, b[i] = a[i  -per/2:i+1+per/2]
	# per is even, left end + [:per/2-1], right end + [-per/2:]
	# per is  odd, left end + [:per/2  ], right end + [-per/2:]
	array = np.array(array)
	atype, shape = array.dtype, array.shape
	if (axis < 0) : axis = len(shape) + axis
	if (axis >= len(shape)) : Raise(Exception, 'axis='+str(axis)+', array.shape='+str(shape)+', axis out of array.shape')
	# Move axis to axis=0
	array = ArrayAxis(array, axis, 0, 'move')
	shape0 = array.shape
	# ND to 2D
	if   (len(shape0) == 1) : array = array[:,None]
	elif (len(shape0)  > 2) : array = array.reshape(shape0[0], np.prod(shape0[1:]))
	shape = array.shape
	#--------------------------------------------------

	if (not reduceshape) : 
		if (shape[1] == 1) : Nprocess = 1
		else : 
			Nprocess = NprocessCPU(Nprocess)[0]
			if (Nprocess > shape[1]) : Nprocess = shape[1]
		if (which == 'per') : weight = GaussianFilter(per, times)[0]
		Nw = len(weight)
		#--------------------------------------------------
		nla = Nw/2 if(Nw%2==1)else Nw/2-1
		nra = Nw/2
		aleft  = np.zeros([nla, shape[1]])
		if (nla != 0) : aleft = aleft + array[0:1]
		aright = np.zeros([nra, shape[1]]) + array[shape[0]-1:shape[0]]
		array = np.concatenate([aleft, array, aright], 0)
		aleft = aright = 0 #@
		bcast = [weight, nla, nra, sigma]
		if (Nprocess == 1) : 
			iterable = ((None,None), array.T, bcast)
			array = _Multiprocess_Smooth(iterable)
		else : 
			pool = PoolFor(0, array.shape[1], Nprocess)
			array = pool.map_async(_Multiprocess_Smooth, array.T, bcast)
			array = np.concatenate(array, 1)
		if (sigma is True) : 
			arrstd = array[shape[0]:].reshape(shape0)
			array  = array[:shape[0]]
		array = array.reshape(shape0)
	#--------------------------------------------------

	else : 
		Nprocess = NprocessCPU(Nprocess)[0]
		for i in xrange(times) : 
			shape = array.shape
			if (Nprocess > shape[0]) : Nprocess = shape[0]
			n1 = PoolFor(0, shape[0], shape[0]/per).nsplit
			#--------------------------------------------------
			if (Nprocess == 1) : 
				iterable = (None, [n1, array], sigma)
				array = _Multiprocess_SmoothReduce(iterable)
			#--------------------------------------------------
			else : 
				n2 = PoolFor(0, len(n1), Nprocess).nsplit
				send = []
				for i in xrange(len(n2)) : 
					j, k = n2[i]
					n3 = n1[j:k]
					send.append( [n3, array[n3.min():n3.max()]] )
				pool = PoolFor()
				array = pool.map_async(_Multiprocess_SmoothReduce, send, sigma)
				array = np.concatenate(array, 0)
		#--------------------------------------------------
		if (sigma is True) : 
			arrstd = array[:,shape[1]:]
			array  = array[:shape[1]]
		shape = (array.shape[0],) + shape0[1:]
		array = array.reshape(shape)
	#--------------------------------------------------
	array = ArrayAxis(array, 0, axis, 'move')
	return array










#
#def _Multiprocess_Smooth( iterable ) : 
#	n1, n2 = iterable[0]
#	array = iterable[1].T
#	weight, lw, dnwa, dnwa1, dnwa2 = iterable[2]
#	shapew = np.ones(len(array.shape), int)
#	shapew[0] = lw
#	weight = weight.reshape(shapew)
#	b = array*0.
#	for i in xrange(len(array)) : 
#		if (i < lw/2) : 
#			dn = lw/2-i
#			da = np.concatenate([array[:1] for j in xrange(dn)])
#			ai = np.concatenate([da, array[:lw/2+1+i]])
#		elif (i >= len(array)-lw/2) : 
#			dn = lw/2 - (len(array)-1 - i)
#			da = np.concatenate([array[-1:] for j in xrange(dn)])
#			ai = np.concatenate([array[i-lw/2:], da])
#		else : ai = array[i-lw/2:i+lw/2+1]
#		b[i] = (ai * weight).sum(0)
#	if (dnwa) : b = b[dnwa1:-dnwa2]
#	return b
#
#
#
#def Smooth( array, axis, per, times=1, sigma=False, reduceshape=False, Nprocess=1 ) : 
#	'''
#	Smooth/Average/Mean array along one axis.
#	We can also use spsn.convolve() to do this, but spsn.convolve() will cost much more memory and time, so, the function written here is the best and fastest.
#
#	Weighted average, sigma will reduce to 1/sqrt{per**times}
#	Equivalent to average over per**times
#
#	axis:
#		array will be smoothed/averaged along which axis.
#
#	per:
#		How many bins/elements to average.
#
#	times:
#		How many times to smooth.
#		For random noise, times=4 is OK
#		Note that this parameter is just for reduceshape=False
#
#	sigma:
#		False, True, int, np.array
#		If False, don't return the error.
#		If True, calculate the error from input array.
#		If int or np.array, use this sigma to calculate the error of the result.
#
#	reduceshape:
#		False, return.shape = array.shape
#		True, return.shape < array.shape
#
#	# Also, we can use Convolve to do it:
#	#    w = array*0
#	#    w = w[len(array)/2-per/2:len(array)/2+per/2+1] = 1
#	#    return Convolve(array, w/w.sum())
#
#	Note that:
#		Wide of 2 windows: w1 < w2
#		a1  = Convolve(a,  w1)
#		a2  = Convolve(a,  w2)
#		a12 = Convolve(a1, w2)
#		=> a12 = Smooth(a2)
#		But the beam sizes of the result maps are similar (roughly the same), Beamsize(a12) >= Beamsize(a2).
#	'''
#	if (per%2 == 0) : per +=1
#	if (per<=1 or times<=0) : return array
#	array = np.array(array)
#	per, atype, shape = int(round(per)), array.dtype, array.shape
#	if (axis < 0) : axis = len(shape) + axis
#	if (axis >= len(shape)) : Raise(Exception, 'axis='+str(axis)+', array.shape='+str(shape)+', axis out of array.shape')
#	# Move axis to axis=0
#	array = ArrayAxis(array, axis, 0, 'move')
#	array0 = array.copy()
#	shape0, shape = array.shape, array.shape
#	# ND to 2D
#	if (len(shape0) == 1) : array = array[:,None]
#	elif (len(shape0) > 2) : array = array.reshape(shape0[0], np.prod(shape0[1:]))
#	shape1 = array.shape  # 2D, along axis=0
#	#--------------------------------------------------
#	if (not reduceshape) : 
#		if (shape1[1] == 1) : Nprocess = 1
#		else : 
#			Nprocess = NprocessCPU(Nprocess)[0]
#			if (Nprocess > shape1[1]) : Nprocess = shape1[1]
#		weight = SmoothWindow( per, times )
#		lw = len(weight)
#		#--------------------------------------------------
#		dnwa, dnwa1, dnwa2 = False, None, None
#		if (shape1[0] < lw) : 
#			dnwa = True
#			dnwa1 = (lw-len(array))/2
#			dnwa2 = lw-len(array) - dnwa1
#			if (dnwa1 > 0) : a1 = np.concatenate([array[:1]  for i in xrange(dnwa1)])
#			else : a1 = array[:0]
#			a2 = np.concatenate([array[-1:] for i in xrange(dnwa2)])
#			array = np.concatenate([a1, array, a2])
#			a1 = a2 = 0 #@
#		#--------------------------------------------------
#		if (Nprocess == 1) : 
#			iterable = ((0,0), array.T, [weight, lw, dnwa, dnwa1, dnwa2])
#			b = _Multiprocess_Smooth(iterable)
#		else : 
#			pool = PoolFor(0, shape1[1], Nprocess)
#			b = pool.map_async(_Multiprocess_Smooth, array.T, [weight, lw, dnwa, dnwa1, dnwa2])
#			b = np.concatenate(b, 1)
#		if (len(shape0) != 2) : b = b.reshape(shape0)
#	#--------------------------------------------------
#	else : 
#		if (times == 1) : Npix = len(array)/per
#		else : Npix = times
#		n = np.linspace(0, len(array), Npix+1).astype(int)
#		b = np.zeros((Npix-1,) + shape[1:])
#		for i in xrange(len(n)-1) : 
#			b[i] = array[n[i]:n[i+1]].mean(0)
#	#--------------------------------------------------
#	b = ArrayAxis(b, 0, axis, 'move')
#	return b
#
#
