import os
from npfmt import *
from Path import *



def Pause() : raw_input()



def Purge() : 
	try : os.system('purge')
	except : pass



def Num( a, dtype=None ) : 
	a = npfmt(a, dtype).flatten()[0]
	return a




def CatalogDecSign( catalogname ) : 
	'''
	In astronomy catalog, Dec is from +90 to -90, but for Dec = -00 23 15, if we np.loadtxt() it, it will miss the "-", treat -00 as 00, and make Dec=-0deg23arcmin15arcsec to +0deg23arcmin15arcsec.
	This function is used to get the sign of Dec

	catalogname: 
		the path/name of the input catalog

	Return:
		1D np.array with the same size of the input catalog

	NOTE THAT:
		In this catalog, must only Dec has sign !
	'''
	instr = open( AbsPath(catalogname) ).readlines()
	sign = np.ones(len(instr), np.float32)
	for i in range(len(sign)) : 
		if ('-' in instr[i]) : sign[i] = -1
	return sign



