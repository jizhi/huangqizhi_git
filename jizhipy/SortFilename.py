

def SortFilename( strlist ) : 
	'''
	(1) If filename is started from a number, it will be placed in front of that started from a letter: '21cmmap', 'apple'
	(2) For the number, we take the 'first whole' number, not just the first number: ['a12', 'a0147', 'a6-5'] will be sorted to be ['a6-5', 'a12', 'a0147'], the 'first whole' numbers are 6(not 6-5), 12, 0147=147
	(3) For the letter, sorted by the first letter from 'a' to 'z'
	'''
	# Sort from number to letter
	order = []
	for i in xrange(len(strlist)) : 
		order.append()



def SortFilename( string ) : 
	'''
	a is any type of list or array.
	Take the first 'whole' number, and sort.
	For example:
		a = ['a12', 'b2-3', 'c35', 'd7_0', 'e4']
		SortStrNum(a) = ['b2-3', 'e4', 'd7_0', 'a12', 'c35']
		first 'whole' numbers are 2, 4, 7, 12, 35

	This function is useful for sort filename
	a = [23, 100, 4, 'xy', 'abc0', 'a12', 'b2-3', 'c35', 'd7_0', 'e4']
	'''
	n = []
	astr = np.array(a, str)
	for i in range(len(astr)) : 
		b = astr[i]
		n1, n2 = len(b), len(b)
		for j in range(len(b)) : 
			if (48 <= ord(b[j]) <= 57) : 
				n1 = j
				break
		for j in range(n1+1, len(b)) : 
			if (48 > ord(b[j]) or ord(b[j]) > 57) : 
				n2 = j
				break
		if (n1 == len(b)) : n = n + [-1]
		else : n = n + [int(b[n1:n2])]
	n = np.array(n)
	n[n<0] = n.max()+1
	nmax = n.max()
	n = np.append([n], [np.arange(len(astr))], 0)
	n = Sort(n)
	n1 = -1
	for i in range(len(n[0])) : 
		if (n[0,i] == nmax) : 
			n1 = i
			break
	n = n[1]
	astr = []
	for i in range(len(a)) : 
		astr = astr + [a[n[i]]]
	if (n1 != -1) : 
		astr = astr[:n1] + list(np.sort(astr[n1:]))
	return astr

