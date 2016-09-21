from npfmt import *
from ArrayAxis import *
from ProgressBar import *
from SysMemory import *
from Smooth import *
from Plot import *



def SelectLeastsq( array, axis=0, firstN=None, usememratio=0.5, method='eachother', verbose=False ): 
	'''
	array is N-D array, for example, a.shape=4x5x6
	I want to get a[:,i,:] which is the resonable one along the second axis

	method:
		'std' or 'eachother'
		(1) use 'std' to judge
		(2) use difference of each other to judge

	firstN:
		Return the first N index

	usememratio:
		There are memcan memory can be used
		Here we use usememratio*memcan

	return:
		list of index, array[index] is the least-square
	'''
	if (method.lower() in ['std','mean']) : method = 'std'
	else : method = 'eachother'
	# Move axis to 0
	array = ArrayAxis(npfmt(array), axis, 0, 'move')
	shape = npfmt(array.shape)
	# flatten()
	array = array.flatten().reshape(shape[0], shape[1:].prod())
	norder = np.arange(len(array))
	# firstN
	try : firstN = int(round(firstN))
	except : firstN = 1
	if (firstN > len(norder)) : firstN = len(norder)
	if (firstN == 0) : firstN = 1
	# Compare to the mean
	array -= Smooth(array.mean(0), 0, array.shape[1]/10, 1)
	nmean = (array**2).mean(1)
	nmean = np.sort(nmean + 1j*norder).imag.astype(int)
	N = len(nmean)/2
	if (firstN > N) : N = firstN
	nmean = nmean[:N]
	if (method == 'std') : return nmean[:firstN]
	#----------
	array = array[nmean]
	#----------
	# Memory of array
	memarray = array.size*1e-8*800 # MB
	totmem = SysMemory()
	memcan = (totmem - memarray)*usememratio # MB
	# Size of each block
	BlockSize = int(1e8/800*memcan)
	# Handle number of rows once
	Nonce = BlockSize/array.size  # Once row
	N = np.linspace(0, shape[0], shape[0]/Nonce+1).astype(int)
	if (N.size == 1) : N = np.concatenate([N,[shape[0]]])
	# Result
	which = []
	if (verbose) : progressbar = ProgressBar('SelectLeastsq():', len(N)-1)
	for i in range(len(N)-1) : 
		if (verbose) : progressbar.Progress()
		a = array[N[i]:N[i+1]]*1
		a = (a[:,None] - array[None,:])**2
		shape = npfmt(a.shape)
		a = a.flatten().reshape(shape[0], shape[1:].prod())
		a = a.mean(1)
		which.append(a)
	# Order
	which = np.concatenate(which)
	which = which + 1j*np.arange(which.size)
	which = np.sort(which).imag.astype(int)
	which = which[:firstN]
	''' OLD, slow
	a = ArrayAxis(npfmt(array), axis, 0, 'move')
	b = np.zeros(len(a))
	for i in range(len(a)) : b[i] = ((a-a[i:i+1])**2).mean()
	which = Sort(b + 1j*np.arange(b.size)).imag.astype(int)
	'''
	return nmean[which]
