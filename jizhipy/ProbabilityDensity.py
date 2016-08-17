
from Edge2Center import *



def BinsArcsinh( array, nbins=None ) : 
	'''
	When use ProbabilityDensity() and RemoveBeyond(), we need to set bins.
	This function is used to return a reasonable and un-uniform bins.

	nbins = bins.size
	'''
	# Get the "real" mean of array
	array = npfmt(array)
	if (nbins is None) : nbins = array.size/10
	if (nbins < 5) : nbins = 5
	amin, amax = array.min(), array.max()
	am, sa = array.mean(), 3*array.std()
	am1, am2 = am-sa, am+sa
	a1, a2 = np.sort([am1, am2])
	array = array[(array>a1)*(array<a2)]
	am = array.mean()  #@
	# np.arcsinh bins
	b1, b2 = np.arcsinh(amin-am), np.arcsinh(amax-am)
	bins = np.sinh(np.linspace(b1, b2, nbins)) + am
	bins[0]  = bins[0] -0.01*abs(bins[0])
	bins[-1] = bins[-1]+0.01*abs(bins[-1])
	return bins



def ProbabilityDensity( array, bins=None, density=True ) : 
	'''
	Return the probability density or number counting of array.
	
	array:
		Input array will flatten()
	
	bins:
		bins can be int or list/np.array.
		If bins == int:
			Number of the uniform bins
		If bins == list/np.array:
			"Edges" of the bins. Note that edges.size = 1+bins.size
			You can set any non-uniform bin widths with this parameter, such as log bin and arcsinh bin.
	
	density:
		If True, return the probability density = counting / total number / bin width
		If False, return the counting number of each bin

	Return:
		[x, xc, y]
		x  is the edge of the bins.
		xc is the center of the bins.
		y  is the probability density of each bin, 
	'''
	if (bins is None) : bins = BinsArcsinh(array)
	y, x = np.histogram(array, bins=bins, density=density)
	xc = Edge2Center(x)
	return [x, xc, y]


