import os
from IsType import *



def DirStr( dirstr, expanduser=False ) : 
	'''
	Convert '../paon4/paon4_data' to '../paon4/paon4_data/'
	
	dirstr:
		Can be str or list of str ['', '', ...]

	expanduser:
		True or False: os.path.expanduser('~/...')

	return:
		Same shape as input dirstr
	'''
	istype = IsType()
	if (istype.isstr(dirstr)) : islist, dirstr = False, [dirstr]
	else : islist, dirstr = True, list(dirstr)
	for i in xrange(len(dirstr)) : 
		if (dirstr[i] == '') : 
			if (islist == False) : return ''
		elif (dirstr[i][-1] != '/') : 
			dirstr[i] = str(dirstr[i]+'/')
		else : dirstr[i] = str(dirstr[i])
		if (expanduser): dirstr[i] = os.path.expanduser(dirstr[i])
	if (islist == False) : dirstr = dirstr[0]
	return dirstr

