from IsType import *
from ShellCmd import *



def Escpath( path, which='rm', abspath=False ) : 
	''' 
	Remove/Add escape characters from/to the path

	In Linux:
		in shell and GUI
		(1) Filename can't contain "/" but can contain ":"
		(2) ":" will be shown as ":"
	In Mac OSX:
		in shell: 
			(1) Filename can't contain "/" but can contain ":"
			(2) ":" will be shown as ":", like in Linux
		in GUI:
			(1) Filename can't contain ":" but can contain "/"
			(2) "/" will be shown as "/" in GUI, but will be shown as ":" in shell

	path:
		Relative or absolute paths are OK

	which:
		=='rm'/'remove': remove "\" from the escape character
		=='add': add "\" to the escape character
	'''
	escchar = '()[]{}<>`!@$^&*=|;,? :'
	path = list(path)
	while ('\\' in path) : path.remove('\\')
	if (which.lower() == 'add') : 
		for i in xrange(len(path)) : 
			if (path[i] in escchar or path[i] in ['"', "'"]) : 
				path[i] = '\\' + path[i]
	if (path == []) : path = ''
	else : 
		for i in xrange(1, len(path)) : path[0] += path[i]
		path = path[0]
		if (path[-1] == '/') : path = path[:-1]
	if (abspath) : path = os.path.abspath(os.path.expanduser(path))
	return path



##################################################
##################################################
##################################################



def Finddir( dirpath, exclude=[], name='' ) : 
	''' Shell command: find dirpath -name name, plus+ exclude '''
	istype = IsType()
	if (not istype.isstr(name)) : name = ''
	name = Escpath(name, 'add')
	#--------------------------------------------------
	dirpath = Escpath(dirpath, 'add', True)
	if (name) : walk = ShellCmd('find '+dirpath+' -name '+name)
	else : walk = ShellCmd('find '+dirpath)
	#--------------------------------------------------
	exctmp = []
	if (istype.isstr(exclude)) : exclude = [exclude]
	for i in xrange(len(exclude)) : 
		exclude[i] = Escpath(exclude[i], 'add')
		if ('*' not in exclude[i]) : 
			dirname=os.path.abspath(os.path.expanduser(exclude[i]))
			exctmp += ShellCmd('find '+dirname)
		else : 
			dirname, basename = os.path.split(exclude[i])
			dirname = os.path.abspath(os.path.expanduser(dirname))
			exctmp += ShellCmd('find '+dirname+' -name '+basename)
	exclude = exctmp
	#--------------------------------------------------
	for i in xrange(len(exclude)) : 
		try : walk.remove(exclude[i])
		except : pass
	walk.remove(dirpath)
	dirs, files = [], []
	for i in xrange(len(walk)) : 
		if (os.path.isdir(walk[i])) : dirs.append(walk[i][len(dirpath)+1:])
		else : files.append(walk[i][len(dirpath)+1:])
	#--------------------------------------------------
	return [dirpath, dirs, files]



