
import numpy as np



def Jy2K( S, DorResol, freq=None, ksd=None, kind=None ) : 
	'''
	Convert flux density S(Jy) to brightness temperature T(K) with Gaussian Beam with FWHM resolution or Rayleigh_criterion resolution (base on the DorResol).

	S:
		the flux density in Jy at this frequency (S will change with frequency). 1Jy = 1e-26 w/m^2/Hz

	DorResol:
		Diameter or resolution in rad.

	freq:
		If DorResol is set to be Diameter, you don't need to set the freq.
		But if DorResol is set to be resolution, you must also set the freq in MHz.

	ksd:
		Beam of antenna maybe not a perfect Gaussian beam.
		sigma = lambda/D /ksd, ksd=1.962~2.286
		Default, ksd=None, 1.03
		For real case, ksd=2

	kind: 
		Just for DorResolu is Diameter.
		We will use D to calculate the solid angle. But for the solid angle, use "FWHM" or "Rayleigh_criterion" ?

	return T in K.

	* S: actually is Jy/beam, because here, we assume the size of the image of the source is equal to the FHWM of Gaussian beam.

	* \Omega = 1.1244*FWHM^2 is very fit by NVSS paper!! Not the \Omega = pi/4 * FWHM^2
	* The result of ksd=None, kind=None is also fit by that of Planck 2015 paper.

	**** 2D flux to K
	Flux2D at any FWHM.
	beam = BeamModel('gaussian', 'dish', D, freq, fov, N)[0] 
	beam = beam / beam.max()  # center is 1
	Then: 
		T2D = Jy2K(Flux2D, D)
		T2D = Convolve(T2D, beam)  # beam, not beam/beam.sum()
	Or:
		Flux2D = Convolve(Flux2D, beam)
		T2D = Jy2K(Flux2D, D)
		(the result is the same)
	'''
	k = 3.2205e-4
	if (freq is not None) : 
		T = (300./freq)**2 * k * S / DorResol**2
	elif (freq is None) : 
		if ksd is None : ksd = (8*np.log(2))**0.5 /1.03
		#	if (kind=='fwhm' or kind=='FWHM') : k1 = 1.03
		#	elif (kind == 'Rayleigh_criterion') : k1 = 1.22
		#	else : raise Exception(efn()+'name must = "FWHM" or "Rayleigh_criterion"')
	#	if (not(1.962 <= ksd <= 2.287)) : raise Exception(efn()+'ksd='+str(ksd)+' is out of range [1.962 ~ 2.287]')
		T = k * S * DorResol**2 * ksd**2 / (8*np.log(2))
#	else : raise Exception(efn()+'Both kind="" and freq="" is not allowed')
	return T
		
