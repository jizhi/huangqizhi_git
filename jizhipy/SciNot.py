from npfmt import *
from IsType import *



def Print( array, precision=6, suppress=True ) : 
	'''
	Format the printing of np.array
	suppress:
		=False, print 1.23e+4
		=True,  print 12340.
	'''
	array = np.array(array)
	precision = int(round(precision))
	if (not suppress) : 
		fmt = '%i.%ie' % (8+precision, precision)
		formatter = {'all': lambda x: format(x, fmt)}
	#	formatter = {'all': lambda x: (('%'+fmt) % x)}
	else : formatter = None
	default = np.get_printoptions()
	np.set_printoptions(precision=precision, suppress=suppress, formatter=formatter)
	print array
	np.set_printoptions(**default)



##################################################
##################################################
##################################################



def SciNot( array ) : 
	'''
	Scientific Notation.
	value can be scale(int/float), list/n-D array
	Return [a, n], value = a * 10**n
	'''
	istype = IsType()
	if (istype.isint(array) or istype.isfloat(array)) : islist = False
	else : islist = True
	array = npfmt(array)
	# Flatten
	shape = array.shape
	array = array.flatten()
	# Get the sign
	sign = np.sign(array)  # sign(0)=0
	# Convert to abs
	array = abs(array)
	# because sign(0)=0, convert 0 to 1
	array[array==0] = 1
	nlarge, nsmall = (array>=1), (array<1)  # bool, not int
	# Use log10 to get the power index
	# >=1
	if (nlarge.sum() > 0) : idxlarge = np.log10(array[nlarge]).astype(int)
	else : idxlarge = []
	# <1
	if (nsmall.sum() > 0) : 
		scalesmall = int(round(np.log10(array[nsmall].min())))-2
		array[nsmall] /= 10.**scalesmall
		idxsmall = np.log10(array[nsmall]).astype(int) + scalesmall
		array[nsmall] *= 10.**scalesmall
	else : idxsmall = []
	# valid and idx
	idx = np.zeros(array.size, int)
	idx[nlarge], idx[nsmall] = idxlarge, idxsmall
	valid = sign * (array / 10.**idx)
	valid, idx = valid.reshape(shape), idx.reshape(shape)
	if (islist) : return (valid, idx)
	else : return (valid[0], idx[0])
	
	
