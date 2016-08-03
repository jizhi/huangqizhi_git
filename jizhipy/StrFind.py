from IsType import *



def StrFind( string, find ) : 
	'''
	find, string:
		Both str. Where is (str)find in the (str)string? Return the index range

	find:
		=Any str that will be found
		='NumType', where is the number(int or float)?

	For example:
		StrFind('CygA665S1dec15', 'S') => [(7,8)]
		It means string[7:8] == 'S'

	return: list of tuple, even just one element
	'''
	#--------------------------------------------------
	istype = IsType()
	if (find != 'NumType') : 
		nfind = []
		if (not istype.isstr(find)) : find = str(find)
		for i in xrange(len(string)-len(find)+1) : 
			if (string[i:i+len(find)] == find) : 
				nfind.append((i,i+len(find)))
		return nfind
	#--------------------------------------------------
	else : 
		strnum = ['' for i in xrange(len(string))]
		nnum, n, nfind, strfind = strnum[:], 0, [], []
		for i in xrange(len(string)) : 
			if (string[i] in '0123456789') : 
				strnum[n] += string[i]
				nnum[n] += str(i)+','
			elif (string[i] == '.') : 
				try : 
					if ((string[i-1] in '0123456789') and (string[i+1] in '0123456789')) : 
						strnum[n] += string[i]
						nnum[n] += str(i)+','
				except : pass
			else : n +=1
		for i in xrange(len(strnum)) : 
			if (strnum[i] != '') : 
				strfind.append(strnum[i])
				nnumi = nnum[i].split(',')
				nfind.append((int(nnumi[0]), int(nnumi[-2])+1)) 
		return [nfind, strfind]


