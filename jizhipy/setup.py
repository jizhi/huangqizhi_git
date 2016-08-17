#! /usr/bin/env python
import os
import sys
from ShellCmd import *


opt = sys.argv[1:]

if (len(opt) == 0) : 
	print 'Miss options, do nothing'
	exit()


syspath = sys.path[:]
site = ''
for i in xrange(len(syspath)) : 
	if (syspath[i][-13:]=='site-packages' and site=='') : site = syspath[i]+'/'
	syspath[i] += '/'


install = False
usrpath, bashrc = '', False
for i in xrange(len(opt)) : 
	if (opt[i].lower() == 'install') : install = True
	else : 
		usrpath = os.path.abspath(os.path.expanduser(opt[i]))
		if (usrpath[-1] != '/') : usrpath += '/'
		if (usrpath not in syspath) : bashrc = True
if (usrpath == '') : usrpath = site


if (install) : 
	try : 
		os.mkdir(usrpath)
		os.system('cp -r ../testpy '+usrpath+'testpy')
	except : 
		os.system('sudo mkdir '+usrpath)
		os.system('sudo cp -r ../testpy '+usrpath+'testpy')


if (bashrc) : 
	pythonpath = 'PYTHONPATH="'+usrpath
	uname = ShellCmd('uname -a')[0][:5]
	if   (uname == 'Darwi') : bashrc = '~/.bash_profile'
	elif (uname == 'Linux') : bashrc = '~/.bashrc'
	bashrc = os.path.expanduser(bashrc)
	if (not os.path.exists(bashrc)) : 
		txt = open(bashrc, 'w')
		print >> txt, pythonpath+'"\n'
	else : 
		tf = False
		txt = open(bashrc, 'r').readlines()
		for i in xrange(len(txt)) : 
			if (txt[i][:10] == 'PYTHONPATH') : 
				tf = True
				break
		if (tf) : pythonpath = pythonpath + ':$PYTHONPATH'
		txt.append( '\n'+pythonpath+'"\n' )
		open(bashrc, 'w').writelines(txt)







		
		


