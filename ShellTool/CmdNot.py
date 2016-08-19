#! /usr/bin/env python
'''
Function:
	Provide 2 commands: 
		cpnot: cp a directory except some files
		rmnot: rm a directory except some files

Usage, examples:
	cpnot -not "a b c" X Y
	rmnot -not "a b c" X Y
'''


import sys
import os


a = sys.argv[1:]

operate = a[0]
a = a[1:]

exclude = a[1+a.index('-not')]
a.remove('-not')
a.remove(exclude)
try : a.remove('-r')
except : pass

exclude = exclude.split(' ')
while ('' in exclude) : exclude.remove('')

b = '-name '+exclude[0]
for i in xrange(1, len(exclude)) : 
	b += ' -or -name '+exclude[i]

if (not os.path.exists(a[1])) : os.mkdir(a[1])

cmd = 'find '+a[0]+' -mindepth 1 -not \( '+b+' \) | xargs -I files '+operate+' -r files '+a[1]
os.system(cmd)
