from PoolFor import *
from ArrayAxis import *
from Same import *
from Raise import *
from Same import *



def _DoMultiprocess_ResetMasked( iterable ) : 
	n1, n2 = iterable[0][0]
	mask, which = iterable[2]
	if (which) : nless, less = iterable[1]
	else : nmore, more = iterable[1]
	nrange = np.arange(mask.size)
	nmasknidx, nmasknvalue = [], []
	for i in xrange(n2-n1) : 
		if (which) : 
			masktmp = mask[nless[i]+2:mask.size]
			j = np.where(masktmp==False)[0][0] + nless[i]+2
			nmasknidx.append(less[i])
			less[i] = j-nless[i]-1
		else : 
			masktmp = mask[:nmore[i]-2+1][::-1]
			j = nmore[i]-2 - np.where(masktmp==False)[0][0]
			nmasknidx.append(more[i])
			more[i] = nmore[i]-j-1
		nmasknvalue.append(nrange[j])
	lm = less if(which)else more
	return [npfmt(nmasknidx), npfmt(nmasknvalue), lm]



def _DoMultiprocess_ResetMasked_npix( iterable ) : 
	n1, n2 = iterable[0][0]
	npix, data = iterable[1]
	value = []
	for i in range(n2-n1) : 
		value.append(np.zeros(npix[i])+data[i])
	return np.concatenate(value)





def ResetMasked( maskedarray, axis, Nprocess=None ) : 
	'''
	Use the elements around the masked elements to calculate some values, then use these values to fill/reset the masked elements
	Return a non-masked array

	maskedarray:
		Any dimension np.ma.MaskedArray()

	axis:
		Along which axis to guess the values

	return:
		resetvalue = ResetMasked(...)
		Usage: maskedarray[maskedarray.mask] = resetvalue
	'''
	if (maskedarray.mask.sum() == 0) : return np.array([])
	Nprocess = NprocessCPU(Nprocess, verbose=False)[0]
	mask0, mask = maskedarray.mask.copy(), maskedarray.mask
	#--------------------------------------------------
	if (axis == 0) : maskedarray = maskedarray.T
	else : 
		maskedarray = ArrayAxis(maskedarray.data, axis, -1, 'move')
		mask = ArrayAxis(mask, axis, -1, 'move')
		maskedarray = np.ma.MaskedArray(maskedarray, mask)
	shape = npfmt(maskedarray.shape)
	maskedarray = maskedarray.flatten()
	#--------------------------------------------------
	nperiod = np.arange(maskedarray.size) % shape[-1]
	nmaskn = nperiod[maskedarray.mask]-1 # left of mask in period
	nperiod = 0 #@
	nrange = np.arange(maskedarray.size)
	#--------------------------------------------------
	less = np.arange(nmaskn.size)[nmaskn<0] # less than 0, other row, the less the better
	nmaskn = nrange[maskedarray.mask]-1 # convert from period to normal order
	nless = nmaskn[less]
	#--------------------------------------------------
	if (Nprocess <= 1) : 
		for i in xrange(len(nless)) : 
			masktmp=maskedarray.mask[nless[i]+2:maskedarray.size]
			j = np.where(masktmp==False)[0][0] + nless[i]+2
			nmaskn[less[i]] = nrange[j]
			less[i] = j-nless[i]-1
		masktmp = 0 #@
	#--------------------------------------------------
	elif (len(nless) > 0) : 
		pool = PoolFor(0, len(nless), Nprocess)
		retn = pool.map_async(_DoMultiprocess_ResetMasked, send=(nless, less), bcast=(maskedarray.mask, True))
		nmasknidx, nmasknvalue, less = [], [], []
		for i in xrange(len(retn)) : 
			nmasknidx.append(retn[i][0])
			nmasknvalue.append(retn[i][1])
			less.append(retn[i][2])
		nmasknidx = np.concatenate(nmasknidx)
		nmasknvalue = np.concatenate(nmasknvalue)
		less = np.concatenate(less)
		nmaskn[nmasknidx] = nmasknvalue
		del pool, retn, nmasknidx, nmasknvalue
	#--------------------------------------------------
	nless = nless+1 + 1j*less
	maskn = maskedarray.mask.copy()
	maskn[nmaskn] = True
	maskn[maskedarray.mask] = False
	nmaskntmp = nrange[maskn]
	#--------------------------------------------------
	# Repeated
	nrepeat = Same(nmaskn)
	for i in range(len(nrepeat)) : 
		n = np.where(nmaskntmp==nrepeat[i,0])[0][0]
		nmaskntmp = np.concatenate([nmaskntmp[:n], np.ones(nrepeat[i,1]-1,int)*nrepeat[i,0], nmaskntmp[n:]])
	nmaskn = nmaskntmp
	#--------------------------------------------------
	nperiod = nrange % shape[-1]
	nmaskp = nperiod[maskedarray.mask]+1
	nperiod = 0 #@
	more = np.arange(nmaskp.size)[nmaskp>=shape[-1]]
	nmaskp = nrange[maskedarray.mask]+1
	nmore = nmaskp[more]
	#--------------------------------------------------
	if (Nprocess <= 1) : # nmaskp, more
		for i in range(len(nmore)) : 
			masktmp = maskedarray.mask[:nmore[i]-2+1][::-1]
			j = nmore[i]-2 - np.where(masktmp==False)[0][0]
			nmaskp[more[i]] = nrange[j]
			more[i] = nmore[i]-j-1
		masktmp = 0 #@
	#--------------------------------------------------
	elif (len(nmore) > 0) : 
		pool = PoolFor(0, len(nmore), Nprocess)
		retn = pool.map_async(_DoMultiprocess_ResetMasked, send=(nmore, more), bcast=(maskedarray.mask, False))
		nmasknidx, nmasknvalue, more = [], [], []
		for i in xrange(len(retn)) : 
			nmasknidx.append(retn[i][0])
			nmasknvalue.append(retn[i][1])
			more.append(retn[i][2])
		nmasknidx = np.concatenate(nmasknidx)
		nmasknvalue = np.concatenate(nmasknvalue)
		more = np.concatenate(more)
		nmaskp[nmasknidx] = nmasknvalue
		del pool, retn, nmasknidx, nmasknvalue
	#--------------------------------------------------
	nmore = nmore-1 + 1j*more
	maskp = maskedarray.mask.copy()
	maskp[nmaskp] = True
	maskp[maskedarray.mask] = False
	nmaskptmp = nrange[maskp]
	#--------------------------------------------------
	# Repeated
	nrepeat = Same(nmaskp)
	for i in range(len(nrepeat)) : 
		n = np.where(nmaskptmp==nrepeat[i,0])[0][0]
		nmaskptmp = np.concatenate([nmaskptmp[:n], np.ones(nrepeat[i,1]-1,int)*nrepeat[i,0], nmaskptmp[n:]])
	nmaskp = nmaskptmp
	nrange = 0 #@
	#--------------------------------------------------
	npix = nmaskp - nmaskn -1
	npix[npix==-1] = np.sort(np.concatenate([nless,nmore])).imag.astype(int)
	#--------------------------------------------------
	if (Nprocess <= 1) : 
		value = []
		for i in xrange(len(npix)) : 
			value.append(np.zeros(npix[i])+(maskedarray.data[nmaskn[i]]+maskedarray.data[nmaskp[i]])/2.)
	elif (len(npix) > 0) : 
		data = (maskedarray.data[nmaskn] + maskedarray.data[nmaskp]) /2.
		pool = PoolFor(0, len(npix), Nprocess)
		value = pool.map_async(_DoMultiprocess_ResetMasked_npix, send=(npix, data))
		del pool
	value = np.concatenate(value)
	npix = nmaskn = nmaskp = 0 #@
	#--------------------------------------------------
	maskedarray.data[maskedarray.mask] = value
	maskedarray = maskedarray.data.reshape(shape)
	if (axis == 0) : maskedarray = maskedarray.T
	else : maskedarray = ArrayAxis(maskedarray, -1, axis, 'move')
	return maskedarray[mask0]


