#! /usr/bin/env python
import os
import sys
import time


os.system('git status')
print '--------------------------------------------------'
print
print 'Do you really want to add and commit all of these?'
print 'Press Y to yes, others to no.'

a = raw_input()

if (a.lower() in ['y', 'yes']) : 
	os.system('git add .')
	mdtime=time.strftime('%Y/%m/%d %p %I:%M:%S', time.localtime())
	os.system('git commit -m "'+mdtime+'"')
	os.system('git push -u origin master')




