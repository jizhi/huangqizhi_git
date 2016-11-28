from Raise import *



def OrderKwargs( stack=2 ) : 
	code = SysFrame(stack, stack)[-1][0]
	n1, n2 = code.find('('), code.rfind(')')
	code = code[n1+1:n2].split(',')
	for i in xrange(len(code)) : 
		c = code[i]
		c = c[:c.find('=')]
		if (c == '') : code[i] = c
		else : 
			for j in xrange(len(c)) : 
				if (c[j] != ' ') : break
			for k in xrange(len(c)-1, 0, -1) : 
				if (c[k] != ' ') : break
			code[i] = c[j:k+1]
	return code

