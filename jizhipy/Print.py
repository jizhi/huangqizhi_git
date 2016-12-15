import sys
import warnings



def Print( a, precision=6, suppress=True ) : 
	'''
	Format the printing of np.array
	suppress:
		=False, print 1.23e+4
		=True,  print 12340.
	'''
	import numpy as np
	np.set_printoptions(precision=precision, suppress=suppress)
	print a





stdoutback = sys.stdout
stderrback = sys.stderr
def PrintFilter( tf=False ) : 
	if (tf) : 
		class Std( object ) : 
			def write( self, stream ) : pass
		sys.stdout = Std()
		sys.stderr = Std()
	else : 
		sys.stdout = stdoutback
		sys.stderr = stderrback





warningfilter, numwarning = [], []
def WarningFilter( tf=None ) : 
	'''
	Print Warning or not?
	tf:
		(1) =None: reset to default
		(2) =True: ignore / not print warning
		(3) =False: print warning
	'''
	if (numwarning == []) : 
		warningfilter.append( warnings.filters[:] )
	numwarning.append(1)
	#--------------------------------------------------
	exceptionwarning = False
	always, ignore = [], []
	for i in xrange(len(warnings.filters)) : 
		w = warnings.filters[i]
		if ('warning' in w[2].__name__.lower()) : 
			if   (w[0] == 'always') : always.append(i)
			elif (w[0] == 'ignore') : 
				ignore.append(i)
				if (w[2].__name__ == 'Warning') : 
					exceptionwarning = True
	#--------------------------------------------------
	if (tf is None) : 
		warnings.filters = warningfilter[0][:]
		return exceptionwarning
	#--------------------------------------------------
	elif (tf) : 
		for i in xrange(len(always)) : 
			warnings.filters[always[i]] = ('ignore',) + warnings.filters[always[i]][1:]
		warnings.filterwarnings('ignore')
	#--------------------------------------------------
	else : 
		for i in xrange(len(ignore)) : 
			warnings.filters[ignore[i]] = ('always',) + warnings.filters[ignore[i]][1:]
	return tf




