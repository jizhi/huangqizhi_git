#! /usr/bin/env python
import os
import sys
from ShellCmd import *
##################################################


# dest, without / at the end
whichpython = ShellCmd('which python')[0].split('/')
whichpython = '/'+whichpython[1]+'/'

syspath = sys.path[:]
for i in xrange(len(syspath)) : 
	if (syspath[i][:len(whichpython)]==whichpython and syspath[i][-9:]=='-packages') : 
		syspath = syspath[i] + '/'
		break

syspath += 'jizhipy'
usrpath = '~/.jizhipy'

whoami = ShellCmd('whoami')[0]
if (whoami == 'root') : 
	dest = syspath
	envlist = []
else : 
	dest = usrpath
	# system environment
	envlist = ['export PYTHONPATH="'+os.path.abspath(os.path.expanduser(usrpath))+'"']


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
# rm setup.py, ShellCmd.py, *.pyc
files.remove('setup.py')
pwd = ShellCmd('pwd')[0].split('/')[-1]
if (pwd == 'ShellTool') : files.remove('ShellCmd.py')
elif ('__init__.py' not in files): os.system('touch __init__.py')
n = 0
for i in xrange(len(files)) : 
	if (files[i-n][-4:] == '.pyc') : 
		files.pop(i-n)
		n +=1


# install / uninstall
if (dest[-1] == '/') : dest = dest[:-1]
if (which == 'uninstall') : 
	print 'Uninstall from  '+dest
	opt = '-rm'
	dest = os.path.abspath(os.path.expanduser(dest+'/'))
	os.system('rm -rf '+dest+'*')
#	for i in xrange(len(files)) : 
#		f = dest + files[i]
#		if (os.path.exists(f)) : os.system('rm -rf '+f)
#		if (os.path.exists(f+'c')) : os.system('rm -rf '+f+'c')

elif (which == 'install') : 
	print 'Install to  '+dest
	opt = '-add'
	dest = os.path.abspath(os.path.expanduser(dest+'/'))
	if (not os.path.exists(dest)) : os.mkdir(dest)
	else : os.system('rm -rf '+dest+'*')
	for i in xrange(len(files)) : 
		os.system('cp -r '+files[i]+' '+dest)


# Add/Remove environment
if (type(envlist) == str) : envlist = [envlist]
for i in xrange(len(envlist)) : 
	os.system('sysenv '+opt+' '+envlist[i])

