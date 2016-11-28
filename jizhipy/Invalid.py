import numpy as np



def Invalid( array, which=True ) : 
	'''
	nan, inf are invalid value.

	which
		True: 
			just mask them, shape of array won't change
			return MaskedArray
		False: 
			return the valid value, will flatten the array
			return normal array
		An value OR ndarray: 
			the masked will be set to be this value
			return MaskedArray
	'''
	array = np.ma.masked_invalid(array)
	if (which is True) : pass
	elif (which is False) : array = array.data[True-array.mask]
	else : array.data[array.mask] = which
	return array


