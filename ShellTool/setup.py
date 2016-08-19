#! /usr/bin/env python
import os
import sys
from ShellCmd import *
##################################################


# dest
whoami = ShellCmd('whoami')[0]
if (whoami == 'root') : dest = '/usr/bin'
else : dest = '~/bin'


# system environment
envlist = ['alias cpnot="python '+dest+'/CmdNot.py'+' cp"', 
'alias rmnot="python '+dest+'/CmdNot.py'+' rm"', 
'alias sysenv="python '+dest+'/SysEnvironment.py"']


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
if (which == 'uninstall') : 
	print 'Uninstall from  '+dest
	opt = '-rm'
	dest = os.path.expanduser(dest+'/')
	for i in xrange(len(files)) : 
		f = dest + files[i]
		if (os.path.exists(f)) : os.system('rm -r '+f)
		if (os.path.exists(f+'c')) : os.system('rm -r '+f+'c')

elif (which == 'install') : 
	print 'Install to  '+dest
	opt = '-add'
	dest = os.path.expanduser(dest+'/')
	if (not os.path.exists(dest)) : os.mkdir(dest)
	for i in xrange(len(files)) : 
		os.system('cp -r '+files[i]+' '+dest)


# Add/Remove environment
if (os.path.exists('/usr/bin/SysEnvironment.py')) : sysenv = '/usr/bin/SysEnvironment.py'
elif (os.path.exists(os.path.expanduser('~/bin/SysEnvironment.py'))) : sysenv = os.path.expanduser('~/bin/SysEnvironment.py')
elif (os.path.exists('SysEnvironment.py')) : sysenv = 'SysEnvironment.py'

if (type(envlist) == str) : envlist = [envlist]
for i in xrange(len(envlist)) : 
	os.system('python '+sysenv+' '+opt+' '+envlist[i])
