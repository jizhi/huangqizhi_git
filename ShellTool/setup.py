#! /usr/bin/env python
import os
import sys

def ShellCmd( cmd ) : 
	strlist = os.popen(cmd).readlines()
	for i in xrange(len(strlist)): strlist[i] = strlist[i][:-1]
	return strlist

##################################################


# dest
whoami = ShellCmd('whoami')[0]
if (whoami == 'root') : dest = '/usr/bin'
else : dest = '~/bin'


##################################################


try : which = sys.argv[1].lower()
except : 
	print 'Error: must be "setup.py install" or "setup.py uninstall"'
	exit()
if (which not in ['install', 'uninstall']) : 
	print 'Error: must be "setup.py install" or "setup.py uninstall"'
	exit()


# Files to be handled
files = ShellCmd('ls')
n = 0
for i in xrange(len(files)) : 
	if ('.py' in files[i-n][-4:]) : 
		files.pop(i-n)
		n +=1


# install / uninstall
if (dest[-1] == '/') : dest = dest[:-1]

if (which == 'uninstall') : 
	print 'Uninstall from  '+dest
	dest = os.path.abspath(os.path.expanduser(dest+'/'))
	for i in xrange(len(files)) : 
		f = dest + files[i]
		if (os.path.exists(f)) : os.system('rm -r '+f)

elif (which == 'install') : 
	print 'Install to  '+dest
	dest = os.path.abspath(os.path.expanduser(dest+'/'))
	if (not os.path.exists(dest)) : os.mkdir(dest)
	for i in xrange(len(files)) : 
		os.system('cp -r '+files[i]+' '+dest)


