from npfmt import *
##################################################


class NoiseSource( object ) : 


	def __init__( self, pixstart=None, pixlength=None, pixperiod=None ) : 
		''' NOTE THAT here set the pixel/index, not second ! 
		Generately, we set pixstart of the first hdf5(hdf5list[0]). However, because inttime is not an integer, after long time, we may stager on pixel, in this case, we set pixstart of current (or one before current) Antarray (see self.Mask() below), then we may have a better result
		'''
		self.pixstart, self.pixlength, self.pixperiod = int(round(pixstart)), int(round(pixlength)), int(round(pixperiod))


	def Mask( self, antarray ) : 
		''' antarray: class:AntArray '''
		Nt = antarray.vis.shape[0]
		nhdf5 = antarray.Hdf5.nhdf5
	#	Nhdf5 = len(antarray.Hdf5.hdf5list)
		Nhdf5 = nhdf5 + 2
		#------------------------------
		# Check pixstart
		if (self.pixstart >= (nhdf5+1)*Nt) : self.pixstart -= Nt
		#------------------------------
		nstart = np.arange(self.pixstart, Nhdf5*Nt, self.pixperiod)
		nend   = np.arange(self.pixstart+self.pixlength, Nhdf5*Nt+1, self.pixperiod)
		if (nend[-1]<nstart[-1]): nend = np.append(nend, [Nhdf5*Nt])
		nstart = nstart[(nhdf5*Nt<=nstart)*(nstart<(nhdf5+1)*Nt)]
		nend   = nend[(nhdf5*Nt<=nend)*(nend<=(nhdf5+1)*Nt)]
		if (nstart[0]>nend[0]) : nstart = np.append(nstart, [nstart[0]-self.pixlength])
		if (nstart[-1]>nend[-1]) : nend = np.append(nend, [nend[-1]+self.pixlength])
		if (nstart[0]<nhdf5*Nt) : nstart[0] = nhdf5*Nt
		if (nend[-1]>(nhdf5+1)*Nt) : nend[0] = nhdf5*Nt
		nstart, nend = nstart-nhdf5*Nt, nend-nhdf5*Nt
		#------------------------------
		self.nmask = np.append(nstart[:,None], nend[:,None], 1)
		# Each nmask[i], vis[nmask[i,0]:nmask[i,1]] is the range of noise source
		self.mask = np.zeros([Nt,], bool)
		for i in xrange(len(self.nmask)) : 
			self.mask[self.nmask[i,0]:self.nmask[i,1]] = True
		# self.mask.shape = (vis.shape[0],)
		# len(self.nmask) = Total number of Noise source
		# These are for current antarray
		return self.nmask, self.mask


##################################################
##################################################
##################################################


class Masking( object ) : 


	def __init__( self, antarray=None ) : 
		if (antarray is None) : return
		self.antarray = antarray  #@ will it increast the memory?

		# Always have self.mask.shape == vis.shape
		shape = list(self.antarray.vis.shape)
		shape[2] = len(self.antarray.visorder)
		self.mask = np.zeros(shape, bool)


	def MaskNoiseSource( self, pixstart, pixlength, pixperiod ) :
		'''
		self.noisesource.pixstart
		self.noisesource.pixlength
		self.noisesource.pixperiod
		self.noisesource.nmask
		self.noisesource.mask
		self.masknoisesource
		'''
		self.noisesource = NoiseSource(pixstart, pixlength, pixperiod)
		self.noisesource.Mask(self.antarray)
		self.masknoisesource = self.noisesource.mask[:,None,None]
		self.mask += self.masknoisesource


	def MaskLoop( self, timeper=60, timetimes=1, freqper=3, freqtimes=1, nsigma=2, nloop=3 ) : 
		'''
		timeper, timetimes, freqper, freqtimes:
			Use for smooth()

		nsigma:
			Calculate sigma of each freq and vis along time (sigma.shape=(Nf, Nv)), value > nsigma*sigma will be considered as RFI

		nloop:
			Each time we will mask the value>nsigma*sigma, and then recalculate a new sigma, and do again. nloop set how many times we will do.

		return self.maskloop
		Total mask: self.maskloop + self.masknoisesource
		'''
		# Read the whole real data as MaskedArray
		vis = np.ma.MaskArray(self.antarray.vis[:,:,self.antarray.visorder], self.mask)
		# Convenient to Smooth(), using non-masked array is much faster than masked array
		vis.data[vis.mask] = vis.mean() # data
		while (nloop > 0) : 
			# time
			dvis = abs(vis).data - Smooth(abs(vis).data, 0, timeper, timetimes) # 3600x512x528
			sigma = dvis.std(0) # 512x528
			vis.mask += (dvis > nsigma*sigma) # mask
			# freq
			dvis = abs(vis).data - Smooth(abs(vis).data, 1, freqper, freqtimes)
			sigma = dvis.std(1)
			vis.mask += (dvis > nsigma*sigma[:,None,:]) # mask
			nloop -=1
		try : self.maskloop = vis.mask - self.masknoisesource
		except : self.maskloop = vis.mask
		self.mask += vis.mask


	def MaskManual( self, maskmanual ) : 
		if (maskmanual.shape != self.mask.shape) : Raise(Exception, 'maskmanual.shape != self.mask.shape')
		self.maskmanual = maskmanual
		self.mask += maskmanual


##################################################
##################################################
##################################################

