from Edge2Center import *
from IsType import *
from Gaussian import *





def BinsNonuniform( array, nbins, root=None, nsigma=None ) : 
	# nsigma
	array = npfmt(array)
	mean, std = array.mean(), array.std()
	if (nsigma is not None) : 
		array = array[(mean-nsigma*sigma<=array)*(array<=mean+nsigma*sigma)]
	mean, std, amin, amax = array.mean(), array.std(), array.min(), array.max()
	# bins
	b = Edge2Center(np.linspace(amin, amax, nbins+1))
	b = GaussianValue(b, mean, std)
	b = np.log(b/b.min()*2)
	try : 
		if (root<=0 or root>1) : root = 1
		b = b**root
	except : pass
	b = 1/b
	b /= b.sum()
	b = np.append([0], np.cumsum(b))
	bins = amin + (amax-amin) * b
	return bins





def ProbabilityDensity( array, bins, nonuniform=False, nsigma=4, density=True ) :
	'''
	Return the probability density or number counting of array.
	
	array:
		Input array will flatten()

	bins, arcsinh:
		(1) bins = list/ndarray (edges), use this, nonuniform is invalid
		(2) bins = int_number
				nonuniform = False : uniform bins
				else : use BinsNonuniform(array, bins, nonuniform), nonuniform is root
	
	nsigma:
		Use how much sigma of the array to calculate the pdf?
		Throw away the points very far from the mean

	density:
		If True, return the probability density = counting / total number / bin width
		If False, return the counting number of each bin

	Return:
		[xe, xc, y]
		xe is the edge of the bins.
		xc is the center of the bins.
		y  is the probability density of each bin, 
	'''
	# Throw away the points very far from the mean
	try : float(nsigma)
	except : nsigma = 4
	array = npfmt(array).flatten()
	sigma, mean = array.std(), array.mean()
	array = array[(mean-nsigma*sigma<=array)*(array<=mean+nsigma*sigma)]
	istype = IsType()
	if (istype.isint(bins) or istype.isfloat(bins)) : 
		bins = int(round(bins))
		if (nonuniform) : bins = BinsNonuniform(array, bins, nonuniform)
		else : bins = np.linspace(array.min(), array.max(), bins+1)
	y, x = np.histogram(array, bins=bins, density=density)
	xc = Edge2Center(x)
	return [x, xc, y]


