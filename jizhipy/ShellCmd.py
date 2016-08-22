from DirStr import *



def mkdir( path=None ) : 
	'''
	If path doesn't exist, mkdir it
	path: str of list of str

	For example:
		mkdir('ab/cd/ef/gh/ij/')
		will do: mkdir('ab/'), mkdir('ab/cd/'), mkdir('ab/cd/ef/'), mkdir('ab/cd/ef/gh/'), mkdir('ab/cd/ef/gh/ij/')
	'''
	istype = IsType()
	if (path is None) : return
	if (istype.isstr(path)) : path, islist = [path], False
	else : path, islist = list(path), True
	path = DirStr(path, True)
	for i in xrange(len(path)) : 
		if (path[i] in ['', './', '../']) : pass
		else : 
			psplit = path[i][:-1].split('/')
			pfont = ''
			for j in xrange(len(psplit)) : 
				pfont = pfont + psplit[j] + '/'
				if (not os.path.exists(pfont)) : os.mkdir(pfont)
	if (not islist) : path = path[0]
	return path



def ShellCmd( cmd ) : 
	'''
	Run shell command and return the result/value

	return:
		strlist
	'''
	if (cmd[:6] == 'mkdir ') : return mkdir(cmd[6:])
	elif (cmd[:3] == 'cd ') : 
		cmd = cmd[3:].split(' ')
		for i in xrange(len(cmd)) : 
			if (cmd[i] != '') : break
		os.chdir(cmd[i])
		return
	strlist = os.popen(cmd).readlines()
	for i in xrange(len(strlist)): strlist[i] = strlist[i][:-1]
	return strlist
	

