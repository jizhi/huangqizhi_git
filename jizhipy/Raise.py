import sys



def SysFrame( downstack=0, upstack=None ) : 
	'''
	downstack: begining of the stack. Default=0 (local)
	upstack: the uppest stack you want. Default to the uppest
	Tack from downstack[include] to upstack[exclude]
	'''
	# How deep of the stack
	for d in xrange(1000000) : 
		try : sys._getframe(d)
		except : break
	stackmax, stackmin = d-1, 0
	if (downstack is None) : downstack = stackmin
	elif (downstack < stackmin) : downstack = stackmin
	if (upstack is None) : upstack = stackmax
	elif (upstack > stackmax) : upstack = stackmax
	if (downstack == upstack) : upstack += 1
	#--------------------------------------------------
	outstr, files, funcs = '', [], []
	for i in range(upstack, downstack, -1) : 
		f = sys._getframe(i)
		filename = f.f_code.co_filename
		lineno   = f.f_lineno
		name     = f.f_code.co_name
		outstr += '  File "'+filename+'", line '+str(lineno)+', in '+name+' =>\n'
		files.append(filename)
		if (name[0] == '<') : funcs.append('')
		else : funcs.append(name)
	outstr = outstr[:-4]
	return [outstr, files, funcs]



def Raise( which=None, message='' ) : 
	'''
	Usage:
		(1) Warning, but run continue
		(2) Exception, stop at once
	'''
	if (which is None) : which = 'Exception'
	if (which in [Warning,'Warning','warning','WARNING']) : 
		which = 'Warning'
		print '-------------------- Warning --------------------'
	else : 
		which = 'Exception'
		print '------------------- Exception -------------------'
	print SysFrame(1)[0]
	print which+': '+str(message)
	print '-------------------------------------------------'
	if (which == 'Exception') : exit()


