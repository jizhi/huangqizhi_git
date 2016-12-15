from IsType import *
from ShellCmd import *





def EscPath( path, which='rm', abspath=False ) : 
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





def FindDir( dirpath, exclude=[], name='' ) : 
	''' Shell command: find dirpath -name name, plus+ exclude '''
	istype = IsType()
	if (not istype.isstr(name)) : name = ''
	name = EscPath(name, 'add')
	#--------------------------------------------------
	dirpath = EscPath(dirpath, 'add', True)
	if (name) : walk = ShellCmd('find '+dirpath+' -name '+name)
	else : walk = ShellCmd('find '+dirpath)
	#--------------------------------------------------
	exctmp = []
	if (istype.isstr(exclude)) : exclude = [exclude]
	for i in xrange(len(exclude)) : 
		exclude[i] = EscPath(exclude[i], 'add')
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





def AbsPath( path ) : 
	'''
	Return the absolute path
	If is directory, '/' at the end
	Else if is file, no '/'
	'''
	path = os.path.abspath(os.path.expanduser(path))
	if (os.path.isdir(path)) : path += '/'
	return path





def ExistsPath( path, old=False ) : 
	'''
	return : 
		If exists, return True, otherwise, return False

	old:
		True: rename the file with '_old'
		False | None: don't rename
	'''
	tf = os.path.exists(path)
	if (old) : 
		n = path.rfind('.')
		if (n < 0) : pathold = path + '_old'
		else : pathold = path[:n] + '_old' + path[n:]
		if (tf) : os.system('mv '+path+' '+pathold)
	return tf


