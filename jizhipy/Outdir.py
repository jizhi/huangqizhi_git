import os
from Raise import *
from ShellCmd import *



def Outdir( stack, name ) : 
	'''
	stack:
		From 0 to uppest, local stack=0
		If you want the uppest, set stack very large (1000000)
		

	name:
		'file' or 'func' ('func' include function and class)

	Outdir(0, 'file') : current file
	Outdir(1, 'file') : 1-upper file
	Outdir(2, 'file') : 2-upper file

	Outdir(0, 'func') : current function/class
	Outdir(1, 'func') : 1-upper function/class
	Outdir(2, 'func') : 2-upper function/class

	return:
		outdir with '/' at the end
	'''
	stack += 1
	# frame
	Nmax = len(SysFrame()[1])-1
	if (stack > Nmax) : stack = Nmax
	frame = SysFrame(stack, stack)
	if (name.lower() in ['func', 'function']) : 
		name = frame[2][0]
		if (name == '') : name = frame[1][0]
	else : name = frame[1][0]
	# name
	n = name.rfind('.')
	if (name[n:] in ['.py', '.pyc']) : name = name[:n]
	if (name[:2] == './') : name = name[2:]
	pwd1 = os.path.abspath('.')
	pwd2 = os.path.abspath('..')
	pwd3 = os.path.abspath('~')
	if (name[:len(pwd1)] == pwd1) : 
		outdir = name[len(pwd1)+1:]
	elif (name[:len(pwd2)] == pwd2) : 
		outdir = '../'+name[len(pwd1)+1:]
	elif (name[:len(pwd3)] == pwd3) : 
		outdir = '~/'+name[len(pwd1)+1:]
	else : outdir = name
	outdir += '_output/'
	outdirabs = os.path.abspath(os.path.expanduser(outdir))
	if (not os.path.exists(outdirabs)) : os.mkdir(outdirabs)
#	return [outdir, outdirabs+'/']
	return outdir


