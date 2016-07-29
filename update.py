#! /usr/bin/env python
import os
import sys
import time


def ShellCmd( cmd ) : 
	strlist = os.popen(cmd).readlines()
	for i in xrange(len(strlist)): strlist[i] = strlist[i][:-1]
	return strlist


# deleted
status = ShellCmd('git status')
delete = []
for i in xrange(len(status)) : 
	if (status[i][1:9] == 'deleted:') : 
		delete = status[i][9:].split(' ')[-1]
		os.system('git rm '+delete)
	

os.system('git status')
print '-------------------------------------------------------'
print
print 'Do you really want to add and commit all of these?'
print 'Press Y to yes, others to no.'

a = raw_input()

if (a.lower() in ['y', 'yes']) : 
	os.system('git add .')
	mdtime=time.strftime('%Y/%m/%d %p %I:%M:%S', time.localtime())
	os.system('git commit -m "'+mdtime+'"')
	os.system('git push -u origin master')
else : print 'Cancel'




