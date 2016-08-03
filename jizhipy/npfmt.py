import numpy as np



def npfmt( value, dtype=None ) : 
	''' convert any value to numpy.array([]) format '''
	if (type(value) == np.ma.core.MaskedConstant) : 
		value = np.ma.asarray([0.])
		value.mask = [True]
	elif (type(value) == np.ma.core.MaskedArray) : pass
	else : 
		try : 
			value = np.array(value, dtype)
			if (value.shape == ()) : value = np.array([value])
		except : value = np.array(value, np.object)
	return value

