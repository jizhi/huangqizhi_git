import jizhipy as jp
from jizhipy.Plot import *
from AntArray import *
from Masking import *
##################################################


class CaliTemperature( object ) : 


	def __init__( self, caligain=None, caliphase=None ) : 
		


	def __init__( self, freq, antarraylist=None, maskinglist=None, caligainlist=None, caliphaselist=None ) :
		'''
		freq:
			in MHz, which freq do you select		
			Must be scale, not list

		antarraylist:
			If we want to calibrate 24 hours data which has 24 hdf5 filts (1 hour per file), we need to mask all files, calibrate the gain and the phase of all files
			antarraylist is [antarray1, antarray2, antarray3, ...]

		maskinglist, caligainlist:
			The same as antarraylist

		caliphaselist:
			(1) One caliphase: all antarraylist use the same caliphase
			(2) List of caliphase [caliphase1, caliphase2, ...], must len(caliphase)==len(antarraylist)
		'''
		istype = jp.IsType()
		n1 = n2 = n3 = n4 = 0
		#--------------------------------------------------
		if (antarraylist is not None) : 
			if (not istype.islist(antarraylist) and not istype.istuple(antarraylist)) : antarraylist = [antarraylist]
			self.antarray = antarray
		try : n1 = len(self.antarray) # Must exist antarraylist
		except : jp.Raise(Exception, 'Missing  self.antarraylist')
		#--------------------------------------------------
		if (maskinglist is not None) : 
			if (not istype.islist(maskinglist) and not istype.istuple(maskinglist)) : maskinglist = [maskinglist]
			self.maskinglist = maskinglist
			n2 = len(self.maskinglist)
			if (n2 != n1) : jp.Raise(Exception, 'len: antarraylist='+str(n1)+' != maskinglist='+str(n2))
		#--------------------------------------------------
		if (caligainlist is not None) : 
			if (not istype.islist(caligainlist) and not istype.istuple(caligainlist)) : caligainlist = [caligainlist]
			self.caligainlist = caligainlist
			n3 = len(self.caligainlist)
			if (n3 != n1) : jp.Raise(Exception, 'len: antarraylist='+str(n1)+' != caligainlist='+str(n3))
		#--------------------------------------------------
		if (caliphaselist is not None) : 
			if (not istype.islist(caliphaselist) and not istype.istuple(caliphaselist)) : caliphaselist = [caliphaselist]
			self.caliphaselist = caliphaselist
			n4 = len(caliphaselist)
			if (n4 == 1) : n4 = len(self.antarray)
			if (n4 != n1) : jp.Raise(Exception, 'len: antarraylist='+str(n1)+' != caliphaselist='+str(n4))
		#--------------------------------------------------
		if (freq is not None) : 
			self.freq = np.array(freq).take(0)
			nfreq = abs(self.antarraylist[0].Ant.freq - self.freq)
			self.nfreq = np.where(nfreq==nfreq.min())[0][0]
		else : jp.Raise(Exception, 'Must set "freq"')
		#--------------------------------------------------
		self.outdir = jp.Outdir((None,'file'), (0,'file'))



	def Calibrator( self, fluxdensity=None, sourcename=None ) : 
		'''
		Just need fluxdensity or sourcename (one of two)

		fluxdensity:
			in Jy, at self.freq
			Use this first

		sourcename:
			Format of the sourcename, CygA as an example:
				sourcename = 'CygA' or 'Cygnus_A'
		'''
		if (fluxdensity is not None) : self.flux = fluxdensity
		elif (sourcename is not None) : 
			bs = jp.BrightSource()
			self.flux = bs.FluxDensity(sourcename, self.freq)



	def Vis( self ) : 
		'''
		Merge antarraylist
		Take the visibility
		nauto = self.antarray.Blorder.auto.size
		self.vis[:,:nauto] is auto
		self.vis[:,nauto:] is cross1
		'''
		vis = []
		for i in xrange(len(self.antarraylist)) : 
			visi = self.antarraylist[i].vis[:,self.nfreq,self.antarraylist[i].visorder]
			if (len(visi.shape) == 1) : visi = visi[:,None]
			vis.append( visi )
		vis = np.concatenate(vis, 0)
		#--------------------------------------------------
		if ('caligainlist' in self.__dict__.keys()) : 
			gaint = []
			for i in xrange(len(self.antarraylist)) : 
				gaint.append(self.caligainlist[i].gaint[:,self.nfreq])
			gaint = np.concatenate(gaint, 0)
		#--------------------------------------------------
			




		antarray = AntArray(self.antarray.Hdf5.hdf5dir)
		auto, cross1 = [], []
		for i in xrange(self.nhdf5list.size) : 
			antarray.WhichHdf5(self.nhdf5list[i])
			auto.append( antarray.vis[:,self.nfreq,self.antarray.Blorder.auto] )
			cross1.append( antarray.vis[:,self.nfreq,self.antarray.Blorder.cross1] )
		auto, cross1 =np.concatenate(auto), np.concatenate(cross1)
		if (len(cross1.shape)==1) : cross1 = cross1[:,None]
		self.vis = np.concatenate([auto, cross1], 1)
		self.vis = np.ma.MaskedArray(self.vis, np.zeros(self.vis.shape, bool))
		self.visorder = np.concatenate([self.antarray.Blorder.auto, self.antarray.Blorder.cross1])
		self.antarray.visorder = self.visorder
		self.masking = Masking(self.antarray)
		self.masking.mask = self.masking.mask[:,0] # (3600,Nv)
		self.masking.mask = np.concatenate([self.masking.mask for i in xrange(self.nhdf5list.size)])
		self.masking.__dict__.pop('antarray')
		#--------------------------------------------------
		self.strbl = list(npfmt(self.antarray.Blorder.blorder[self.antarray.Blorder.auto][:,0], str))
		strc = npfmt(self.antarray.Blorder.blorder[self.antarray.Blorder.cross1], str)
		for i in xrange(len(strc)) : 
			self.strbl.append(strc[i,0]+'-'+strc[i,1])
		#--------------------------------------------------
		Nt = self.antarray.vis.shape[0]
		inttime = self.antarray.Ant.inttime
		n0 = self.nhdf5list[0]*Nt
		self.timem =np.arange(n0, n0+self.vis.shape[0])*inttime/60.
		np.save(self.outdir+'timem.npy', self.timem)



	def MaskNoiseSource( self, pixstart, pixlength, pixperiod ) : 
		if ('vis' not in self.__dict__.keys()) : self.Vis()
		antarray = AntArray(self.antarray.Hdf5.hdf5dir)
		antarray.visorder = self.visorder
		mask = []
		for i in xrange(self.nhdf5list.size) : 
			antarray.WhichHdf5(self.nhdf5list[i])
			masking = Masking(antarray)
			masking.mask = masking.mask[:,0:1,:] # 3D
			masking.MaskNoiseSource(pixstart, pixlength, pixperiod)
			mask.append(masking.masknoisesource.copy()[:,0])
		self.masking.masknoisesource=np.concatenate(mask)
		self.masking.mask += self.masking.masknoisesource
		self.vis.mask += self.masking.mask


	def MaskLoop( self, timeper=60, timetimes=1, freqper=3, freqtimes=1, nsigma=2, nloop=3 ) : 
		if ('vis' not in self.__dict__.keys()) : self.Vis()
		antarray = AntArray(self.antarray.Hdf5.hdf5dir)
		antarray.visorder = self.visorder
		mask = []
		for i in xrange(self.nhdf5list.size) : 
			antarray.WhichHdf5(self.nhdf5list[i])
			antarray.vis = antarray.vis[:,self.nfreq:self.nfreq+1]
			masking = Masking(antarray)
			masking.mask = self.masking.mask[:,None,:] # 3D
			masking.MaskLoop(timeper,timetimes, freqper,freqtimes)
		mask = np.concatenate(mask)[:,0]
		try : self.masking.maskloop = mask - self.masking.masknoisesource
		except : self.masking.maskloop = mask
		self.masking.mask += self.masking.maskloop
		self.vis.mask += self.masking.mask


	def MaskManual( self, maskmanual ) : 
		if ('vis' not in self.__dict__.keys()) : self.Vis()
		if (maskmanual.shape != self.vis.shape) : Raise(Exception, 'maskmanual.shape != self.vis.shape')
		self.masking.maskmanual = npfmt(maskmanual, bool)
		self.masking.mask += self.masking.maskmanual
		self.vis.mask += self.masking.mask


	def Gaint( self, timeper=60, timetimes=10, legendsize=None, markersize=5 ) : 
		''' Actually the sigma of noise of visibility changing with time, but can represent the profile of G(t)

		perpix:
			Every how many time-pixel to calculate one std
		'''
		class Gaint( object ) : pass
		self.gaint = Gaint()
		vis = ResetMasked(self.vis, 0)
		vis -= Smooth(vis, 0, timeper, timetimes)
		vis = np.ma.MaskedArray(vis, self.vis.mask)
		n = np.linspace(0, len(vis), len(vis)/timeper+1).round().astype(int)
		gaintr, gainti = [], []
		for i in xrange(len(n)-1) : 
			gaintr.append( [vis.real[n[i]:n[i+1]].std(0)] )
			gainti.append( [vis.imag[n[i]:n[i+1]].std(0)] )
		gaintr = np.concatenate(gaintr)
		gainti = np.concatenate(gainti)
		self.gaint.gaint = gaintr + 1j*gainti
		self.gaint.gainpix = Edge2Center(n).round().astype(int)
		#--------------------------------------------------
		nauto = self.antarray.Blorder.auto.size
		self.gaint.gaintmean = self.gaint.gaint[0].real.copy()
		for i in xrange(len(self.gaint.gaintmean)) : 
			a = self.gaint.gaint[:,i].copy()
			if (i < nauto) : a = a.real
			else : a = np.concatenate([a.real, a.imag])
			self.gaint.gaintmean[i] = a[abs(a)>1e-10][2:-2].mean()
		#--------------------------------------------------
		abl = self.antarray.Blorder.blorder[self.antarray.Blorder.auto][:,0]
		cbl = self.antarray.Blorder.blorder[self.antarray.Blorder.cross1]
		self.gaint.sqrt2 = np.zeros(len(cbl))
		for i in xrange(len(self.gaint.sqrt2)) : 
			n1 = np.where(abl==cbl[i,0])[0][0]
			n2 = np.where(abl==cbl[i,1])[0][0]
			a1 = self.gaint.gaintmean[n1]
			a2 = self.gaint.gaintmean[n2]
			c  = self.gaint.gaintmean[i+nauto]
			self.gaint.sqrt2[i] = (a1*a2)**0.5 /c
		#--------------------------------------------------
		# Plot gaintmean, sqrt2
		plt.figure(figsize=(17,6))
		plt.subplot(1,2,1)
		color = plt_color(len(self.strbl))
		for i in xrange(len(self.strbl)) : 
			plt.plot(i+1, self.gaint.gaintmean[i], marker='o', ls='', markersize=markersize, color=color[i], label=self.strbl[i])
		plt.legend(fontsize=legendsize)
		plt.xlim(0, len(self.strbl)+1)
		plt.xticks([0], [''])
		plt.ylabel(r'$\sigma$ [A.U.]', size=20)
		plt.title(r'Average $\sigma$ of noise', size=20)
		plt.subplot(1,2,2)
		ncross = len(self.visorder) - nauto
		color = plt_color(ncross)
		plt.plot(np.arange(ncross+2), np.zeros(ncross+2)+2, 'k--')
		for i in xrange(ncross) : 
			plt.plot(i+1, self.gaint.sqrt2[i]**2, marker='o', ls='', markersize=markersize, color=color[i], label=self.strbl[i+nauto])
		plt.legend(fontsize=legendsize)
		plt.xlim(0, ncross+1)
		plt.xticks([0], [''])
		plt.ylim(1.7, 2.3)
		plt_axes('y', 'both', [0.1,0.01], '%.1f')
		plt.ylabel('Ratio', size=20)
		plt.title(r'Ratio = $\sigma_i \cdot \sigma_j \,/\, \sigma_{ij}^2$', size=20)
		plt.savefig(self.outdir+'gaint.png')
		plt.close()


	def CaliTemp( self, sourceflux, kc=0.85 ) : 
		Ts = Jy2K(sourceflux, self.antarray.Ant.dishdiam*0.9)
		nauto = self.antarray.Blorder.auto.size
		self.vis.data[:,:nauto] /= self.gaint.gaintmean[:nauto] # stdauto = 1
	#	self.vis.data[:,nauto:] /= (self.gaint.gaintmean[nauto:]*2**0.5) # stdcross = 1/2**0.5
		self.vis.data[:,nauto:] /= (self.gaint.gaintmean[nauto:]*self.gaint.sqrt2)
		As = abs(self.vis[:,nauto:]).max(0).mean()
		self.calitemp = np.zeros(self.vis.shape[1])
		self.calitemp[:nauto], self.calitemp[nauto:] = (Ts/As *kc), (Ts/As)
		self.Tsys = Ts/As
		self.vis.data[:] *= self.calitemp  # K
		np.save(self.outdir+'calitemp', self.calitemp)
		np.save(self.outdir+'vis_K.npy', self.vis.data)
		np.save(self.outdir+'vis_K.mask.npy', self.vis.mask, bool)


	def Smooth( self, times=0 ) : 
		vis = ResetMasked(self.vis, 0)  # return real array
		vis = Smooth(vis, 0, 3, times)
		self.vis = np.ma.MaskedArray(vis, self.vis.mask)
		self._smoothdone = True


	def Plot( self, plotnv=None, pltaxes=None ) : 
		# Plot all auto
		nauto = self.antarray.Blorder.auto.size
		ymin, ymax = abs(self.vis[:,:nauto]).min(), abs(self.vis[:,:nauto]).max()
		dy = ymax-ymin
		ymin = ymin-0.1*dy if(ymin>0.1*dy)else 0
		ymax += 0.1*dy
		plt.figure(figsize=(17,6))
		plt.subplot(1,2,1)
		color = plt_color(nauto)
		for i in xrange(nauto) : 
			if (self._smoothdone) : 
				plt.plot(self.timem, self.vis.data[:,i].real, color=color[i], label=self.strbl[i])
			else : 
				plt.plot(self.timem, self.vis[:,i].real, color=color[i], label=self.strbl[i])
		plt.xlim(self.timem.min(), self.timem.max())
		plt.ylim(ymin, ymax)
		if (pltaxes is not None) : pltaxes
		plt.xlabel(r't [min]', size=16)
		plt.ylabel(r'T [K]', size=16)
		plt.title('Auto-correlations', size=16)
		#--------------------------------------------------
		plt.subplot(1,2,2)
		if (plotnv is None) : # longest East-West baseline
			bl = self.antarray.Blorder.baseline[self.antarray.Blorder.cross1]
			bl = abs(bl[:,0])
			nv = np.where(bl==bl.max())[0][0]
		else : nv = plotnv
		if (self._smoothdone) : 
			plt.plot(self.timem, self.vis.data[:,nv+nauto].real, 'b-', label=self.strbl[nv+nauto]+', real')
			plt.plot(self.timem, self.vis.data[:,nv+nauto].imag, 'r-', label=self.strbl[nv+nauto]+', imag')
		else : 
			plt.plot(self.timem, self.vis[:,nv+nauto].real, 'b-', label=self.strbl[nv+nauto]+', real')
			plt.plot(self.timem, self.vis[:,nv+nauto].imag, 'r-', label=self.strbl[nv+nauto]+', imag')
		plt.xlim(self.timem.min(), self.timem.max())
		plt.ylim(-1.05*dy, 1.05*dy)
		if (pltaxes is not None) : pltaxes
		plt.xlabel(r't [min]', size=16)
		plt.ylabel(r'T [K]', size=16)
		plt.title('Cross-correlation of baseline='+self.strbl[nv+nauto], size=16)
	#	plt.suptitle(r'Calibrat the amplitude to Kelvin @ '+str(self.freq)+'MHz, $T_{sys}='+str(self.Tsys)+'$ K', size=20)
		plt.suptitle(r'Calibrat the amplitude to Kelvin @ '+str(self.freq)+'MHz', size=20)
		plt.savefig(self.outdir+'vis_'+self.strbl[nv+nauto]+'_'+str(self.freq)+'MHz.png')
		plt.close()



