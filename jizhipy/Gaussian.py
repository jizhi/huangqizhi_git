from Raise import *
import numpy as np



def GaussianValue( x, mean, std ) : 
	y = 1/(2*np.pi)**0.5 /std *np.e**(-(x-mean)**2 /(2*std**2))
	return y

 
def LogNormalValue( x, mean, std ) : 
	if (x.min() < 0) : Raise(Exception, 'x.min() < 0')
	y = 1/(2*np.pi)**0.5 /std /x *np.e**(-(np.log(x)-mean)**2 /(2*std**2))
	return y

