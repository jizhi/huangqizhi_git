#! /usr/bin/env python
import os
import sys


def ShellCmd( cmd ) : 
	a = os.popen(cmd).readlines()
	for i in xrange(len(a)) : a[i] = a[i][:-1]
	return a


#def Status() : 
#	a = ShellCmd('git status')
#	n1, n2 = None, len(a)-1
#	for i in xrange(len(a)) : 
#		if (a[i][:11] == '  (use "git') : n1 = i+2
#	a = a[n1:n2]
#	for i in xrange(len(a)) : 
#		ai = a[i][1:]
#		if (ai[0] == '.') : continue
		
	

print a
