import numpy as np



def CompactDimension( array ) : 
	array = np.array(array)
	shape = []
	for i in array.shape : 
		if (i != 1) : shape.append(i)
	array = array.reshape(shape)
	return array



