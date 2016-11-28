import os
from npfmt import *



def Pause() : raw_input()



def Purge() : 
	try : os.system('purge')
	except : pass



def Print( a, precision=6, suppress=True ) : 
	'''
	Format the printing of np.array
	suppress:
		=False, print 1.23e+4
		=True,  print 12340.
	'''
	np.set_printoptions(precision=precision, suppress=suppress)
	print a



def Num( a, dtype=None ) : 
	a = npfmt(a, dtype).flatten()[0]
	return a

