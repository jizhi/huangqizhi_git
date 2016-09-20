from npfmt import *
from Raise import *



def ArrayAxis( array, axis1, axis2, act='move' ) : 
	'''
	array:
		Any dimension array (real, complex, masked)

	axis1, axis2:
		move axis1 to axis2

	act:
		'move' or 'exchange'
		move     : move     a[axis1]  to  b[axis2]
		exchange : exchange a[axis1] with a[axis2]
	'''
	array = npfmt(array)
	shapeo = array.shape
	if (len(shapeo) <= 1) : return array
	if (axis1 < 0) : axis1 = len(shapeo) + axis1
	if (axis2 < 0) : axis2 = len(shapeo) + axis2
	if (axis1 == axis2) : return array
	if (axis1>=len(shapeo) or axis2>=len(shapeo)) : Raise(Exception, 'axis1='+str(axis1)+', axis2='+str(axis2)+' out of array.shape='+str(shapeo)+'=>'+str(len(shapeo))+'D')
	if (len(shapeo) == 2) : return array.T
	#--------------------------------------------------
	def _Move( array, axis1, axis2 ) : 
		shapeo = array.shape
		shapen = list(shapeo)[:]
		shapen.remove(shapeo[axis1])
		shapen=tuple(shapen[:axis2]+[shapeo[axis1]]+shapen[axis2:])
		a = np.zeros(shapen, array.dtype)
		for i in xrange(shapeo[axis1]) : 
			if (axis1 == 0) : 
				if   (axis2 == 1) : a[:,i]         = array[i]
				elif (axis2 == 2) : a[:,:,i]       = array[i]
				elif (axis2 == 3) : a[:,:,:,i]     = array[i]
				elif (axis2 == 4) : a[:,:,:,:,i]   = array[i]
				elif (axis2 == 5) : a[:,:,:,:,:,i] = array[i]
			elif (axis1 == 1) : 
				if   (axis2 == 0) : a[i]           = array[:,i]
				elif (axis2 == 2) : a[:,:,i]       = array[:,i]
				elif (axis2 == 3) : a[:,:,:,i]     = array[:,i]
				elif (axis2 == 4) : a[:,:,:,:,i]   = array[:,i]
				elif (axis2 == 5) : a[:,:,:,:,:,i] = array[:,i]
			elif (axis1 == 2) : 
				if   (axis2 == 0) : a[i]           = array[:,:,i]
				elif (axis2 == 1) : a[:,i]         = array[:,:,i]
				elif (axis2 == 3) : a[:,:,:,i]     = array[:,:,i]
				elif (axis2 == 4) : a[:,:,:,:,i]   = array[:,:,i]
				elif (axis2 == 5) : a[:,:,:,:,:,i] = array[:,:,i]
			elif (axis1 == 3) : 
				if   (axis2 == 0) : a[i]           = array[:,:,:,i]
				elif (axis2 == 1) : a[:,i]         = array[:,:,:,i]
				elif (axis2 == 2) : a[:,:,i]       = array[:,:,:,i]
				elif (axis2 == 4) : a[:,:,:,:,i]   = array[:,:,:,i]
				elif (axis2 == 5) : a[:,:,:,:,:,i] = array[:,:,:,i]
			elif (axis1 == 4) : 
				if   (axis2 == 0) : a[i]           = array[:,:,:,:,i]
				elif (axis2 == 1) : a[:,i]         = array[:,:,:,:,i]
				elif (axis2 == 2) : a[:,:,i]       = array[:,:,:,:,i]
				elif (axis2 == 3) : a[:,:,:,i]     = array[:,:,:,:,i]
				elif (axis2 == 5) : a[:,:,:,:,:,i] = array[:,:,:,:,i]
			elif (axis1 == 5) : 
				if   (axis2 == 0) : a[i]         = array[:,:,:,:,:,i]
				elif (axis2 == 1) : a[:,i]       = array[:,:,:,:,:,i]
				elif (axis2 == 2) : a[:,:,i]     = array[:,:,:,:,:,i]
				elif (axis2 == 3) : a[:,:,:,i]   = array[:,:,:,:,:,i]
				elif (axis2 == 4) : a[:,:,:,:,i] = array[:,:,:,:,:,i]
			else : Raise(Exception, 'ArrayAxis() can just handel 6-D array. For >= 7-D array, you can modify this function by yourself.')
		return a
	#--------------------------------------------------
	array = _Move(array, axis1, axis2)
	if (abs(axis1-axis2)!=1 and act.lower()=='exchange') : 
		if (axis1 < axis2) : axis1, axis2 = [axis2-1, axis1]
		else : axis1, axis2 = [axis2+1, axis1]
		array = _Move(array, axis1, axis2)
	return array


