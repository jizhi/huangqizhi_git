import os
import sys



class IsType( object ) : 
	dtype = 'class:IsType'

	def isint( self, a ) : 
		typestr = str(type(a)).split("'")[-2].split('.')[-1][:3]
		if (typestr == 'int') : return True
		else : return False
	
	def isfloat( self, a ) : 
		typestr = str(type(a)).split("'")[-2].split('.')[-1][:5]
		if (typestr == 'float') : return True
		else : return False
	
	def isstr( self, a ) : 
		typestr = str(type(a)).split("'")[-2].split('.')[-1][:3]
		if (typestr == 'str') : return True
		else : return False

	def islist( self, a ) : 
		typestr = str(type(a)).split("'")[-2].split('.')[-1]
		if (typestr == 'list') : return True
		else : return False

	def istuple( self, a ) : 
		typestr = str(type(a)).split("'")[-2].split('.')[-1]
		if (typestr == 'tuple') : return True
		else : return False
	
	def isdict( self, a ) : 
		typestr = str(type(a)).split("'")[-2].split('.')[-1]
		if (typestr == 'dict') : return True
		else : return False
	
	def isnparray( self, a ) : 
		typestr = str(type(a)).split("'")[-2].split('.')[-1]
		if (typestr == 'ndarray') : return True
		else : return False
	
	def isclass( self, a ) : 
		typestr = str(type(a)).split("'")[-2].split('.')[-1]
		if (typestr in ['type', 'classobj']) : return True
		else : return False
		
	def isinstance( self, a ) : 
		typestr1 = str(type(a)).split("'")[-2]
		typestr2 = str(type(a)).split(' ')[0].split('<')[-1]
		if (typestr1=='instance' or typestr2=='class') : 
			return True
		else : return False



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
	

