from npfmt import *
from IsType import *
from Raise import *
from Other import *
from Print import *



def _SymmetryBeam( which, fwhm, theta, normalized ) : 
	'''
	theta, fwhm:
		in rad
		Can be any shape
	
	return:
		beam.shape = (fwhm.shape, theta.shape)

	GaussianBeam = exp{-(theta/FWHM)^2/2/sigma^2}, its max=1, not normalized.

	Set theta=0.5*FWHM, compute exp{-(1/2)^2/2/sigma^2}=0.5(because the max=1, theta=0.5*FWHM will decreate to half power) and get sigma^2=1/(8ln2)

	\int{e^(-a*x^2)dx} = (pi/a)^0.5
	So, for the normalized GaussianBeam_normalized 
	= 1/FWHM * (4ln2/pi)^0.5 * exp(-4ln2 * theta^2 / FWHM^2)

	SincBeam = np.sinc(2.783/np.pi * theta**2 / fwhm**2)
	'''
	if (IsType().isnum(fwhm)) : islistfwhm = False
	else : islistfwhm = True
	if (IsType().isnum(theta)) : islisttheta = False
	else : islisttheta = True
	fwhm, theta = npfmt(fwhm), npfmt(theta)
	shapefwhm, shapetheta = fwhm.shape, theta.shape
	fwhm, theta = fwhm.flatten(), theta.flatten()
	if (which == 'gaussian') : b = np.exp(-4*np.log(2) * theta[None,:]**2 / fwhm[:,None]**2)
	elif (which == 'sinc') : b = np.sinc(2.783/np.pi * theta[None,:]**2 / fwhm[:,None]**2)**2
	if (normalized) : 
		a = 4*np.log(2) / fwhm[:,None]**2
		b = b/(np.pi/a)**0.5
		a = 0 #@
	if (not islistfwhm and not islisttheta) : b = b[0,0]
	elif (not islistfwhm and islisttheta) : 
		b = b.reshape(shapetheta)
	elif (islistfwhm and not islisttheta) : 
		b = b.reshape(shapefwhm)
	else : b = b.reshape(shapefwhm + shapetheta)
	return b



def GaussianBeam( fwhm, theta, normalized=False ) : 
	return _SymmetryBeam('gaussian', fwhm, theta, normalized)


def SincBeam( fwhm, theta, normalized=False ) : 
	return _SymmetryBeam('sinc', fwhm, theta, normalized)





def EllipticGaussianBeam( fwhm1, fwhm2, theta, phi, normalized=False ) : 
	'''
	theta, phi, fwhm1, fwhm2:
		in rad
		theta.shape == phi.shape, can be any shape
		fwhm1.shape == fwhm2.shape, can be any shape
	'''
	if (IsType().isnum(fwhm1)) : islistfwhm1 = False
	else : islistfwhm1 = True
	if (IsType().isnum(fwhm2)) : islistfwhm2 = False
	else : islistfwhm2 = True
	islistfwhm = bool(islistfwhm1 + islistfwhm2)
	if (IsType().isnum(theta)) : islisttheta = False
	else : islisttheta = True
	if (IsType().isnum(phi)) : islistphi = False
	else : islistphi = True
	islisttheta = bool(islisttheta + islistphi)
	fwhm1, fwhm2, theta, phi = npfmt(fwhm1), npfmt(fwhm2), npfmt(theta), npfmt(phi)
	shape1, shape2, shapet, shapep = fwhm1.shape, fwhm2.shape, theta.shape, phi.shape
	printstr = 'fwhm1.shape='+str(shape1)+', fwhm2.shape='+str(shape2)+', theta.shape='+str(shapet)+', phi.shape='+str(shapep)
	if (shape1 != shape2) : Raise(Exception, 'fwhm1.shape != fwhm2.shape. '+printstr)
	if (shapet != shapep) : Raise(Exception, 'theta.shape != phi.shape. '+printstr)
	#--------------------------------------------------
	fwhm1, fwhm2, theta, phi = fwhm1.flatten(), fwhm2.flatten(), theta.flatten(), phi.flatten()
	b = np.exp(-4*np.log(2) * theta[None,:]**2 * ((np.cos(phi[None,:])/fwhm1[:,None])**2 + (np.sin(phi[None,:])/fwhm2[:,None])**2))
	if (normalized) : 
		a = 4*np.log(2) * ((np.cos(phi[None,:])/fwhm1[:,None])**2 + (np.sin(phi[None,:])/fwhm2[:,None])**2)
		b = b/(np.pi/a)**0.5
		a = 0 #@
	#--------------------------------------------------
	if (not islistfwhm and not islisttheta) : b = b[0,0]
	elif (not islistfwhm and islisttheta) : 
		b = b.reshape(shapetheta)
	elif (islistfwhm and not islisttheta) : 
		b = b.reshape(shapefwhm)
	else : b = b.reshape(shapefwhm + shapetheta)
	return b





def ThetaPhiMatrix( thetawidth, Npix ) : 
	'''
	thetawidth:
		total field of view (from left to rignt) in rad.
		Must be one

	Npix:
		theta and phi matrix shape = (Npix, Npix)
		Must be one int

	# At the center (Npix/2, Npix/2), theta = 0

	return:
		[theta, phi] in rad.
	'''
	thetawidth = Num(thetawidth, float)
	Npix = Npix0 = Num(Npix, int)
	if (Npix %2 == 0 ) : Npix += 1
	N = Npix / 2
	# Length of each pixel
	pixlist = np.arange(-N, N+1.)
	dx = thetawidth / (Npix-1)
	theta = dx * (pixlist[:,None]**2 + pixlist[None,:]**2)**0.5
	phi = pixlist[None,:] + 1j*pixlist[::-1,None]
	PrintFilter(True)
	phi = np.arctan(phi.imag / phi.real)
	PrintFilter()
	phi[N,N] = 0
	phi[:,:N] += np.pi
	phi[N+1:,N:] += 2*np.pi
	invalid = (theta > np.pi)
	theta[invalid], phi[invalid] = np.nan, np.nan
	theta, phi = theta[:Npix0,:Npix0], phi[:Npix0,:Npix0]
	return np.array([theta, phi])


