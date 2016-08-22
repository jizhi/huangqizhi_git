import numpy as np



def Invalid( array, mask=True ) : 
	'''
	nan, inf are invalid value.

	mask:
		True    : just mask them, shape of array won't change
		False   : return the valid value, will flatten the array
		An value: the masked will be set to be this value
	'''
	array = np.ma.masked_invalid(array)
	if (mask is True) : pass
	elif (mask is False) : array = array.data[True-array.mask]
	else : 
		array.data[array.mask] = mask
	#	array = array.data
	return array


