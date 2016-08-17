
from npfmt import *
from ShellCmd import *



def GSM(frequency, fmt=None, bit=32, gsmpath='/usr/bin/gsm/'):
	'''
	Generate foreground all sky map by GSM:
		nside=512, ring, K, Galactic.

	frequency: 
		MHz. one value or list or tuple or ndarray

	fmt:
		'npy', 'fits', None
		=None: GSM default save to '.dat'
		='npy': Convert .dat to .npy, then delete .dat file
		='fits': Save to .fits

	bit:
		32, 64
		Save to float32 or float64

	return:
		[[outname, open_function], [,], [,], ...]
		if (type(frequency)=='NumType'): return [outname, open_function]
	'''
	istype = IsType()
	islist = False if(istype.isint(frequency) or istype.isfloat(frequency))else True
	freqlist = npfmt(frequency).flatten()
	outnamelist = []
	openfunclist = []
	for i in freqlist : 
		freqstr = str(i)
		outname = 'gsm_'+freqstr
		outexist = ShellCmd('ls '+outname+'.*')
		if (len(outexist) > 0) : # exist
			if (outname+'.npy' in outexist) : 
				outnamelist.append(outname+'.npy')
				openfunclist.append(np.load)
			elif (outname+'.fits' in outexist) : 
				outnamelist.append(outname+'.fits')
				openfunclist.append(pyfits.getdata)
			else : 
				outnamelist.append(outexist[0])
				openfunclist.append(np.loadtxt)
		else : 
			gsmpath = gsmpath if(gsmpath[-1]=='/')else gsmpath+'/'
			if (not os.path.exists(gsmpath+'gsm_parameter.out')):
				try : 
					gsmpathout = gsmpath
					os.system('gfortran -ffixed-line-length-none '+gsmpath+'gsm_parameter.f -o '+gsmpathout+'gsm_parameter.out')
				except : 
					gsmpathout = ''
					os.system('gfortran -ffixed-line-length-none '+gsmpath+'gsm_parameter.f -o '+gsmpathout+'gsm_parameter.out')
			else : gsmpathout = gsmpath
			os.system(gsmpathout+'gsm_parameter.out '+freqstr+' '+outname+'.dat '+gsmpath)
			if (fmt is None) : fmt = 'dat'
			if (fmt[-3:].lower() == 'npy') : 
				a = np.loadtxt(outname+'.dat')
				if   (bit == 32) : a = np.float32(a)
				elif (bit == 64) : a = np.float64(a)
				np.save(outname, a)
				a = 0 #@
				os.system('rm '+outname+'.dat')
				outnamelist.append(outname+'.npy')
				openfunclist.append(np.load)
			elif (fmt[-4:].lower() == 'fits') : 
				a = np.loadtxt(outname+'.dat')
				if   (bit == 32) : a = np.float32(a)
				elif (bit == 64) : a = np.float64(a)
				Array2FitsImag(a, outname+'.fits')
				a = 0 #@
				os.system('rm '+outname+'.dat')
				outnamelist.append(outname+'.fits')
				openfunclist.append(pyfits.getdata)
			else : 
				outnamelist.append(outname+'.dat')
				openfunclist.append(np.loadtxt)
	returnlist = []
	for i in range(len(outnamelist)) : 
		returnlist.append([outnamelist[i], openfunclist[i]])
	if (not islist) : returnlist = returnlist[0]
	return returnlist


