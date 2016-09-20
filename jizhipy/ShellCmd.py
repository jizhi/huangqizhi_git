import os
import sys
from Raise import *



class cd( object ) : 

	def __init__( self ) : 
		self.originalpath = os.getcwd()

	def go( self, path=None ) : 
		if (path is not None) : os.chdir(path)

	def back( self ) : 
		os.chdir(self.originalpath)



def ShellCmd( cmd ) : 
	'''
	Absolutely use shell command, and return the strings that should be printed on the screen
	return: list of str
	'''
	if (cmd.lower() == 'pdf') : 
		path = SysFrame()[1][-1]
		n = path.rfind('/')
		path = path[:n+1]
		if (os.popen('which xdg-open').readlines()) : which = 'xdg-open'
		if (os.popen('which open').readlines()) : which = 'open'
		os.system(which+' '+path+'python-shellcmd.pdf')
		return
	strlist = os.popen(cmd).readlines()
	for i in xrange(len(strlist)): strlist[i] = strlist[i][:-1]
	return strlist
	

