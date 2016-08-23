import os



def ShellCmd( cmd ) : 
	'''
	Run shell command and return the result/value

	return:
		strlist
	'''
	strlist = os.popen(cmd).readlines()
	for i in xrange(len(strlist)): strlist[i] = strlist[i][:-1]
	return strlist
	

