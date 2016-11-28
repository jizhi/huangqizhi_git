from npfmt import *
from IsType import *



def Freq2Velo( dfreq, freq0 ) : 
	'''
	Doppler effect/shift

	freq0:
		float/on value in MHz
		center frequency

	dfreq:
		in MHz
		can be one/int/float or list/ndarray
		differential from freq0

	return:
		dvelo in km/s (from 0km/s)
		same shape as dfreq
	'''
	freq0 = npfmt(freq0, float).flatten()[0]
	istype, islist = IsType(), True
	if (istype.isint(dfreq) or istype.isfloat(dfreq)) : 
		islist = False
	dfreq = npfmt(dfreq)
	dvelo = 3e5*dfreq/freq0
	if (not islist) : dvelo = dvelo[0]
	return dvelo



def Velo2Freq( dvelo, freq0 ) : 
	'''
	dvelo:
		in km/s
		can be one/int/float or list/ndarray
		differential from 0km/s

	return:
		dfreq in MHz (from freq0)
		same shape as dvelo
	'''
	freq0 = npfmt(freq0, float).flatten()[0]
	istype, islist = IsType(), True
	if (istype.isint(dvelo) or istype.isfloat(dvelo)) : 
		islist = False
	dvelo = npfmt(dvelo)
	dfreq = dvelo*freq0/(3e5)
	if (not islist) : dfreq = dfreq[0]
	return dfreq




