import numpy as np
from IsType import *
from Raise import *



def StrlistAdd( *arg ) : 
	'''
	arg = (str or list of str)
	Combine arg to str or list of str, according to the type of arg (list or not)
	tuple will convert to list

	For example: 
		StrlistAdd('a', 'b') => 'ab'
		StrlistAdd('a', ['b']) => ['ab']
		StrlistAdd('a', ['b'], ['c', 'd']) => ['abc', 'abd']
		StrlistAdd('data/', ['CygA', 'CasA'], '_paon4.fits') => ['data/CygA_paon4.fits', 'data/CasA_paon4.fits']
	'''
	# Check arg
	arg = list(arg)
	narg = np.zeros([len(arg),], int)
	islist = narg*0
	istype = IsType()
	for i in xrange(len(arg)) : 
		if (istype.isstr(arg[i])) : narg[i], islist[i] = 1, 0
		else : narg[i], islist[i] = len(arg[i]), 1
	nmax = narg.max()
	narg = narg[(narg.min()<narg)*(narg<narg.max())]
	if (narg.size > 0) : Raise(Exception, 'Element numbers in arg='+str(arg)+' are not suitable')
	if (islist.max() == 1) : islist = True
	else : islist = False
	#----------
	for i in xrange(len(arg)) : 
		argi = arg[i] if(istype.isstr(arg[i]))else arg[i][0]
		arg[i] = [argi for j in xrange(nmax)]
	restr = []
	for j in xrange(nmax) : 
		strij = ''
		for i in xrange(len(arg)) : strij += arg[i][j]
		restr.append(strij)
	if (islist == False) : restr = restr[0]
	return restr

