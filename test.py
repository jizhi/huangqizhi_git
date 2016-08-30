#! /usr/bin/env python
import jizhipy as jp
import os


def Abspath( path ) : 
	''' 
	Convert to abspath and check whether exists 
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
	'''
	escchar = '()[]{}<>`!@$^&*=\|;,? :'
	which = ''
	for i in xrange(len(escchar)) : 
		if (escchar[i] in path or escchar[i] in ["'",'"']) : 
			which = escchar[i]
			break
	if (which != '') : 
		n = path.find(which)
		pathesc, pathtmp = list(path), ''
		if (path[n-1] != '\\') : 
			for i in xrange(len(path)) : 
				if (pathesc[i] in escchar) : pathesc[i] = '\\' + pathesc[i]
			for i in xrange(len(pathesc)) : pathtmp += pathesc[i]
			pathesc = pathtmp
		else : 
			while ('\\' in pathesc) : pathesc.remove('\\')
			for i in xrange(len(pathesc)) : pathtmp += pathesc[i]
			pathesc = path
			path = pathtmp
	else : pathesc = path
	path = os.path.abspath(os.path.expanduser(path))
	pathesc = os.path.abspath(os.path.expanduser(pathesc))
	exists = os.path.exists(path)
	return [path, pathesc, exists]





a = jp.ShellCmd('find .')
print a
