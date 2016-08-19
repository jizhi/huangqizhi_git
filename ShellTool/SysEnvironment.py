#! /usr/bin/env python
'''
Function:
	Provide 1 command: 
		sysenv: is used to "add" or "remove" the environment in .bashrc or .bash_profile

Usage, examples:
	sysenv -add alias cpnot="~/bin/cmdnot.py cp"
	sysenv -rm  alias cpnot="~/bin/cmdnot.py cp"

	sysenv -add export PYTHONPATH="/python2.7/site-packages/"
	sysenv -rm export PYTHONPATH="/python2.7/site-packages/"
'''


import os
import sys
from ShellCmd import *


# Environment
try : 
	pm = sys.argv[1].lower()
	which = sys.argv[2]
	env = sys.argv[3]
except : 
	print 'Error: SysEnvironment, does nothing'
	exit()
envkey = env[:env.index('=')]
envvalue = env[env.index('=')+1:]
if (envvalue[0] == "'") : envvalue = '"'+envvalue[1:-1]+'"'
if (envvalue[0] != '"') : envvalue = '"'+envvalue+'"'

if (which=='export' and envvalue[-2]=='/') : envvalue = envvalue[:-2]+'"' # rm /


# Which bash
uname = ShellCmd('uname -a')[0][:6]
if (uname == 'Linux ') : bash = '~/.bashrc'
elif (uname == 'Darwin') : 
	if (os.path.exists(os.path.expanduser('~/.bash_profile'))) : bash = '~/.bash_profile'
	elif (os.path.exists(os.path.expanduser('~/.bashrc'))) : bash = '~/.bashrc'
bash = os.path.expanduser(bash)
if (not os.path.exists(bash)) : ShellCmd('touch '+bash)


# List
font = which+' '+envkey+'='
env1 = font+ envvalue[:-1]+'/"' +'\n'  # +" +/
env2 = font+ envvalue +'\n'  # +" -/
env3 = font+ envvalue[1:-1]+'/' +'\n'  # -" +/
env4 = font+ envvalue[1:-1] +'\n'  # -" -/
if (which == 'export') : 
	back = ':$'+envkey
	env5 = font+ envvalue[:-1]+'/'+back+'"' +'\n'  # +" +/
	env6 = font+ envvalue[:-1]+back+'"' +'\n'  # +" -/
	env7 = font+ envvalue[1:-1]+'/'+back +'\n'  # -" +/
	env8 = font+ envvalue[1:-1]+back +'\n'  # -" -/
	envlist = [env1, env2, env3, env4, env5, env6, env7, env8]
else : envlist = [env2, env4]


# Read bash
txt = open(bash).readlines()
while (txt[-1] == '\n') : txt.pop(-1)


# Remove
if (pm.lower() in ['-rm', '-remove', '--rm', '--remove']) : 
	for i in xrange(len(envlist)) : 
		while (envlist[i] in txt) : txt.remove(envlist[i])
	if (which == 'export') : 
		for i in xrange(len(txt)) : 
			if (txt[i][:len(which)+1+len(envkey)] == which+' '+envkey) :
				if (':$'+envkey in txt[i]) : 
					n1 = txt[i].index(':$'+envkey)
					n2 = n1 + len(':$'+envkey)
					txt[i] = txt[i][:n1] + txt[i][n2:]
					break
	while (txt[-1] == '\n') : txt.pop(-1)
	open(bash, 'w').writelines(txt)
	os.system('source '+bash)
	exit()


# Add
if (which == 'export') : envlist = [env1, env5]
envexist = False
for i in xrange(len(envlist)) : 
	if (envlist[i] in txt) : envexist = True
if (envexist) : exit()

if (which == 'alias') : txt.append('\n'+envlist[0])
elif (which == 'export') : 
	keyexist = False
	for i in xrange(len(txt)) : 
		if (txt[i][:len(which)+1+len(envkey)]==which+' '+envkey) :
			keyexist = True
			break
	if (keyexist) : env = env5
	else : env = env1
	txt.append('\n'+env)

open(bash, 'w').writelines(txt)
os.system('source '+bash)
