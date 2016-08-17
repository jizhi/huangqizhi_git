
from ArrayAxis import *



def SmoothWeight( per, times, plv=False ) : 
	if (per%2==0 and plv) : 
		print 'Warning: SmoothWeight(), per must be odd, reset per='+str(per)+' to per='+str(per+1)
		per +=1
	a0 = np.ones(per)
	#----- Method 1 -----
#	for n in xrange(1, times+1) : 
#		if (times == 1) : break
#		if (n == 1) : continue
#		a = np.zeros([per, 1+2*n])
#		for j in xrange(len(a)) : 
#			a[j,j:j+len(a0)] = a0
#		a0 = a.sum(0)
#	weight = a0 / (1.*per)**times
	#----- Method 1 END -----
	#----- Method 2 -----
	for n in xrange(1, times+1) : 
		if (times == 1) : break
		if (n == 1) : continue
		a = np.zeros([per, 1+2*n*(per/2)])
		for j in xrange(len(a)) : 
			a[j,j:j+len(a0)] = a0
		a0 = a.sum(0) /per
	weight = a0/per  #@ Remember this !!!
	#----- Method 2 END -----
	return weight



def Smooth( array, axis, per, times=1, sigma=False, reduceshape=False ) : 
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
		dnwa1 = None
		if (len(array) < lw) : 
			dnwa1 = (lw-len(array))/2
			dnwa2 = lw-len(array) - dnwa1
			if (dnwa1 > 0) : a1 = np.concatenate([array[:1]  for i in xrange(dnwa1)])
			else : a1 = []
			a2 = np.concatenate([array[-1:] for i in xrange(dnwa2)])
			array = np.concatenate([a1, array, a2])
			a1 = a2 = 0 #@
		#--------------------------------------------------
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
		if (dnwa1 is not None) : b = b[dnwa1:-dnwa2]
		ai = w = 0 #@
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




#
#def ArrayLeftRightInterp( array, axis, num, back=False ) : 
#	''' Mostly be used in the function Smooth() '''
#	if (type(array) == np.ma.core.MaskedArray) : return array
#	nalr, array = num, npfmt(array)
#	# move axis to axis=0
#	array = ArrayAxis(array, axis, 0, 'move')
#	if (back == False) : 
#		aleft = array[:nalr]
#		aright = array[-nalr:]
#		aleft = (aleft - aleft[-1:] + aleft[:1])[:-1]
#		array = np.append(aleft, array, axis)
#		aright = (aright - aright[:1] + aright[-1:])[1:]
#		array = np.append(array, aright, axis)
#	else : array = array[nalr-1:-nalr+1]
#	array = ArrayAxis(array, 0, axis, 'move')
#	return array
#
#
#
#def _DoMultiprocess_Smooth( iterable ) : 
#	mp1, mp2 = iterable[0]
#	sigma, reduceshape, per, N, arrayfont, arrayend, array0font, array0end, nsplit, dosigma = iterable[2]
#	#--------------------------------------------------
#	if (sigma is not False) : 
#		arrayloc = iterable[1]
#		arrayloc = ArrayAxis(arrayloc, -1, 0, 'move')
#		nr = len(arrayloc)/2
#		if (dosigma) : arrayloc0 = arrayloc[nr:]
#		arrayloc  = arrayloc[:nr]
#		arrayloc  = ArrayAxis(arrayloc,  0, -1, 'move')
#		if (dosigma) : 
#			arrayloc0 = ArrayAxis(arrayloc0, 0, -1, 'move')
#			sa = arrayloc*0.
#		else : sa = None
#	else : 
#		arrayloc = iterable[1]
#		sa = None
#	#--------------------------------------------------
#	result = arrayloc*0.
#	n = nsplit.index(iterable[0])
#	nf = len(arrayfont[n])
#	ne = nf + len(arrayloc)
#	arrayloc=np.concatenate([arrayfont[n],arrayloc,arrayend[n]])
#	if (sigma is not False and dosigma) : 
#		arrayloc0 = np.concatenate([array0font[n], arrayloc0, array0end[n]])
#	N = len(arrayloc)
#	#--------------------------------------------------
#	for i in xrange(nf, ne) : 
#		if (i < per/2) : n1, n2 = 0, per-per/2+i
#		elif (i > N-(per+1)/2) : n1, n2 = (per+1)/2-per+i, N
#		else : n1, n2 = i-per/2, per-per/2+i
#		#--------------------------------------------------
#		result[i-nf] = arrayloc[n1:n2].mean(0)  # main !
#		#--------------------------------------------------
#		if (sigma is not False and dosigma) : 
#			if (sigma is True) : sa[i-nf] = RMS(arrayloc0[n1:n2]-result[i-nf:i-nf+1], 0) /(n2-n1)
#			else : 
#				if (sigma.size==1) : error = (n2-n1)**0.5*sigma
#				else : error = ((sigma[n1-nf:n2-nf]**2).sum())**0.5
#				sa[i-nf] = error /(n2-n1)
#	#--------------------------------------------------
#	return [result, sa]
#
#
#
#def Smooth( array, axis, per, times=1, sigma=False, reduceshape=False, applr=False, Nprocess=None ) : 
#	'''
#	Smooth/Average/Mean array along one axis.
#	we can also use spsn.convolve() to do this, but spsn.convolve() will spend much more memory and time, so, the function written here is the best and fastest.
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
#	if (per<=1 or times<=0) : return array
#	array = np.array(array)
#	per, atype, shape = int(round(per)), array.dtype, array.shape
#	if (axis < 0) : axis = len(shape) + axis
#	if (axis >= len(shape)) : Raise(Exception, 'axis='+str(axis)+', array.shape='+str(shape)+', axis out of array.shape')
#	# Move axis to axis=0
#	array = ArrayAxis(array, axis, 0, 'move')
#	shape = array.shape  # new shape
#	#--------------------------------------------------
#	if (reduceshape) : applr, times = False, 1
#	#--------------------------------------------------
#	# Append left and right
#	if (applr == True) : 
#		array = ArrayLeftRightInterp(array, 0, per)
#		shape = array.shape
#	#--------------------------------------------------
#	# sigma, True, False, int/ndarray
#	if (sigma is not False) : # True, int/ndarray
#		if (sigma is not True) : 
#			sigma = npfmt(sigma)
#			if (sigma.size == 1) : sigma = array*0.+sigma.take(0)
#			elif (applr == True) : 
#				sigma = ArrayLeftRightInterp(sigma, 0, per)
#		array0 = array*1
#	#--------------------------------------------------
#	#--------------------------------------------------
#	for t in xrange(times) : 
#		if (t == times-1) : dosigma = True
#		else : dosigma = False
#		#--------------------------------------------------
#		pool = PoolFor(0, shape[0], Nprocess)
#		#--------------------------------------------------
#		nsplit = pool.nsplit[:]
#		arrayend, arrayfont = [], []
#		for i in xrange(len(nsplit)) : 
#			arrayfont.append(array[nsplit[i][0]-per:nsplit[i][0]])
#			arrayend.append( array[nsplit[i][1]:nsplit[i][1]+per])
#		if (sigma is False) : array0font, array0end = None, None
#		elif (t == 0) : array0font, array0end = arrayfont[:], arrayend[:]
#		#--------------------------------------------------
#		if (sigma is not False) : # True, int/ndarray
#			if (len(shape) == 1) : array = np.concatenate([array[:,None], array0[:,None]], 1)
#			else : array = np.concatenate([array, array0], -1)
#		#--------------------------------------------------
#		retn = pool.map_async(_DoMultiprocess_Smooth, array, (sigma, reduceshape, per, shape[0], arrayfont, arrayend, array0font, array0end, nsplit, dosigma))
#		#--------------------------------------------------
#		# Distinguish result and sa
#		result, satmp = [], []
#		for i in xrange(len(retn)) : 
#			result.append(retn[i][0])
#			if (retn[i][1] is not None) : satmp.append(retn[i][1])
#		array = np.concatenate(result)
#		if (len(satmp) > 0) : sa = np.concatenate(satmp)
#		del pool, retn, result, satmp
#	#--------------------------------------------------
#	array0 = 0 #@
#	if (applr is True) : 
#		array = ArrayLeftRightInterp(array, 0, per, True)
#		if (sigma is not False) : 
#			sa = ArrayLeftRightInterp(sa, 0, per, True)
#	#--------------------------------------------------
#	if (reduceshape) : 
#		nc = np.linspace(0, shape[0], shape[0]/per+1)
#		nc = Edge2Center(nc).astype(int)
#		array = array[nc]
#		if (sigma is not False) : sa = sa[nc]
#	#--------------------------------------------------
#	array = ArrayAxis(array, 0, axis, 'move')
#	if (sigma is not False) : 
#		sa = ArrayAxis(sa, 0, axis, 'move')
#		return [array, sa]
#	else : return array
#




