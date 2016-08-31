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

usrpath = '~/.python-packages/'
usrpathpwd = os.path.abspath(os.path.expanduser(usrpath))


whoami = ShellCmd('whoami')[0]
if (whoami == 'root') : 
	dest = syspath + 'jizhipy'
	envlist = []
else : 
	dest = usrpath + 'jizhipy'
	# system environment
	home = ShellCmd('echo $HOME')[0]
	if (usrpathpwd[:len(home)] == home) : usrpathpwd = '$HOME'+usrpathpwd[len(home):]
	envlist = ['export PYTHONPATH="'+usrpathpwd+'"']


##################################################


try : which = sys.argv[1].lower()
except : 
	print 'Destination: ', dest
	if (envlist) : print 'Add environment: ', envlist[0]
	print 'setup.py install    to   install'
	print 'setup.py uninstall  to uninstall'
	exit()
if (which not in ['install', 'uninstall']) : 
	print 'Error: must be "setup.py install" or "setup.py uninstall"'
	exit()


# Files to be handled
files = os.listdir('.')  # include .xxx
# rm setup.py, ShellCmd.py, *.pyc
files.remove('setup.py')
if ('__init__.py' not in files): 
	os.system('touch __init__.py')
	files.append('__init__.py')
n = 0
for i in xrange(len(files)) : 
	if (os.path.splitext(files[i-n]) == '.pyc') : 
		files.pop(i-n)
		n +=1


# install / uninstall
if (which == 'uninstall') : 
	print 'Uninstall from  '+dest
	opt = '-rm'
	dest = os.path.abspath(os.path.expanduser(dest))
	os.system('rm -rf '+dest)


elif (which == 'install') : 
	print 'Install to  '+dest
	opt = '-add'
	dest = os.path.abspath(os.path.expanduser(dest))
	if (os.path.exists(dest)) : os.system('rm -rf '+dest)
	os.makedirs(dest)
	for i in xrange(len(files)) : 
		os.system('cp -r '+files[i]+' '+dest+'/')


# Add/Remove environment
for i in xrange(len(envlist)) : 
	env = envlist[i]
	n1 = env.find(' ')
	for n2 in xrange(n1+1, len(env)) : 
		if (env[n2] != ' ') : break
	which = env[:n1]
	env = env[n2:]
	os.system('sysenv '+opt+' '+which+" '"+env+"'")

