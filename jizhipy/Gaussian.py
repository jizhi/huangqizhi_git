from Raise import *
from npfmt import *
from scipy.integrate import quad, dblquad



def GaussianValue( x, mean, std ) : 
	y = 1/(2*np.pi)**0.5 /std *np.exp(-(x-mean)**2 /(2*std**2))
	return y

 

def LogNormalValue( x, mean, std ) : 
	if (x.min() < 0) : Raise(Exception, 'x.min() < 0')
	y = 1/(2*np.pi)**0.5 /std /x *np.exp(-(np.log(x)-mean)**2 /(2*std**2))
	return y



def GaussianSolidAngle( fwhmx, fwhmy=None ) : 
	'''
	fwhmx, fwhmy: 
		in rad

	return: 
		solid_angle = GaussianSolidAngle(fwhmx, fwhmy) in sr
	'''
	fwhmx = npfmt(fwhmx, float).flatten()[0]
	try : fwhmy = npfmt(fwhmy, float).flatten()[0]
	except : fwhmy = fwhmx
	#--------------------------------------------------
	if (fwhmx != fwhmy) : 
		def SolidAngle(y, x) : 
			sa = np.exp(-4*np.log(2) *y**2 *(np.sin(x)**2/fwhmx**2 + np.cos(x)**2/fwhmy**2)) *np.sin(y)
			return sa
		solid_angle = dblquad(SolidAngle, 0., 2*np.pi, lambda x:0., lambda x:np.pi)[0]
	#--------------------------------------------------
	else : 
		def SolidAngle(x) : 
			sa = np.exp(-4*np.log(2) *x**2 /fwhmx**2) *np.sin(x) *2*np.pi
			return sa
		solid_angle = quad(SolidAngle, 0, np.pi)[0]
	#--------------------------------------------------
	return solid_angle


