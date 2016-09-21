from ArrayAxis import *



def Sort( array, along='[0,:]', l2s=False ) : 
	'''
	array:
		Can be any shape

	along:
		Must as format like '[n1,n2,:,n4]'
		Must have ':', use [2,:] instead of [2]
		array[n1,n2,:,n4] must 1D so that we can sort along this
		along=[:,2] : second column

	l2s: 
		l2s=False: from small to large (default)
		l2s=True : from large to small
	'''
	along = along[1:-1].split(',')
	axis = along.index(':')
	along.pop(axis)
	along = np.array(along, int)
	#--------------------------------------------------
	array = npfmt(array)
	if (len(array.shape) == 1) : 
		array = np.sort(array)
		if (l2s) : array = array[::-1]
		return array
	#--------------------------------------------------
	array = ArrayAxis(array, axis, -1, 'move')
	shape = array.shape
	#--------------------------------------------------
	cumprod = np.cumprod((shape[1:-1]+(1,))[::-1])[::-1]
	along = (along*cumprod).sum()
	a = array.reshape(np.prod(shape[:-1]), shape[-1])[along]
	#--------------------------------------------------
	a = a + 1j*np.arange(a.size)
	a = np.sort(a).imag.astype(int)
	if (l2s) : a = a[::-1]
	#--------------------------------------------------
	array = ArrayAxis(array, -1, 0, 'move')
	array = array[a]
	array = ArrayAxis(array, 0, axis, 'move')
	return array


