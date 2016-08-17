from AntArray import *
from Masking import *

import scipy.signal as spsn
from Sph2Circ import *
from Plot import *
from PoolFor import *
from ResetMasked import *
from Smooth import *
from FuncFit import *

##################################################



def _DoMultiprocess_FitBeam( iterable ) : 
	n1, n2  = iterable[0]
	beam = iterable[1]
	x0, t0, s0, Dec, showprogress = iterable[2]
	if (showprogress) : progressbar = ProgressBar('CaliPhase.FitBeam(), pid='+str(os.getpid())+':', len(beam))
	#--------------------------------------------------
	def func(x, p) :  # x is time, not angle
		beta = Circ2Sph((x-p[1])/3600.*15*np.pi/180, Dec*np.pi/180)
		y = p[0] * np.exp(-beta**2 /2 /p[2]**2)
		return y
	#--------------------------------------------------
	pf = []
	for i in xrange(n2-n1) : 
		if (showprogress) : progressbar.Progress()
		p0 = [beam[i].max(), t0, s0]
		pf.append( FuncFit(func, x0, beam[i], p0)[0] )
	return npfmt(pf)


def _DoMultiprocess_FitPhase( iterable ) : 
	n1, n2 = iterable[0]
#	phase = iterable[1]  # smoothdone
#	x0, t0, Dec, freq, bl, showprogress, Nv, lat = iterable[2] # smoothdone
	phase, t0list = iterable[1]
	x0, Dec, freq, bl, showprogress, Nv, lat = iterable[2] # smoothdone
	if (showprogress) : progressbar = ProgressBar('CaliPhase.FitPhase(), pid='+str(os.getpid())+':', len(phase))
	pf = []
	#--------------------------------------------------
	for i in xrange(n2-n1) : 
		t0 = t0list[i] #@
		if (showprogress) : progressbar.Progress()
		nf, nv = (i+n1)/Nv, (i+n1)%Nv
		def func(x, p) : 
			beta =Circ2Sph((x-t0)/3600.*15*np.pi/180, Dec*np.pi/180)
			angle = 2*np.pi/(300./freq[nf]) *( bl[nv,0]*np.sin(beta) - bl[nv,1]*np.cos(beta)*np.sin((Dec-lat)*np.pi/180) ) + p[0]
			return np.sin(angle)
		#--------------------------------------------------
		p0 = np.pi
		pf.append( FuncFit(func, x0, np.sin(phase[i]), p0)[0] )
	return npfmt(pf).flatten()


def _DoMultiprocess_FitVis( iterable ) : 
	n1, n2 = iterable[0]
	vis, A0, t0list = iterable[1]
	x0, s0, Dec, freq, bl, showprogress, Nv, lat = iterable[2]
	if (showprogress) : progressbar = ProgressBar('CaliPhase.FitPhase(), pid='+str(os.getpid())+':', len(phase))
	pf = []
	#--------------------------------------------------
	for i in xrange(n2-n1) : 
		if (showprogress) : progressbar.Progress()
		nf, nv = (i+n1)/Nv, (i+n1)%Nv
		def func(x, p) : 
			beta = Circ2Sph((x-t0list[i])/3600.*15*np.pi/180, Dec*np.pi/180)
			A = p[0] * np.exp(-beta**2 /2 /p[1]**2)
			angle = 2*np.pi/(300./freq[nf]) *( bl[nv,0]*np.sin(beta) - bl[nv,1]*np.cos(beta)*np.sin((Dec-lat)*np.pi/180) ) + p[2]
			return A*np.sin(angle)
		#	return A*np.cos(angle)
		#--------------------------------------------------
		p0 = [A0[i], s0, 1]
		pf.append( FuncFit(func, x0, vis[i], p0)[0] )
	return npfmt(pf)  # shape=(n2-n1, 4)


def _DoMultiprocess_FitExp( iterable ) : 
	n1, n2 = iterable[0]
	vis, t0list = iterable[1]
	x0, Dec, freq, bl, showprogress, Nv, lat = iterable[2]
	if (showprogress) : progressbar = ProgressBar('CaliPhase.FitExp(), pid='+str(os.getpid())+':', len(phase))
	pf = []
	#--------------------------------------------------
	for i in xrange(n2-n1) : 
		if (showprogress) : progressbar.Progress()
		nf, nv = (i+n1)/Nv, (i+n1)%Nv
		beta = Circ2Sph((x0-t0list[i])/3600.*15*np.pi/180, Dec*np.pi/180)
		A = abs(vis[i])
		def func(x, p) : 
			angle = 2*np.pi/(300./freq[nf]) *( bl[nv,0]*np.sin(beta) - bl[nv,1]*np.cos(beta)*np.sin((Dec-lat)*np.pi/180) ) + p[0]
			return A*np.sin(angle)
		#	return A*np.cos(angle)
		#--------------------------------------------------
		p0 = np.pi
		pf.append( FuncFit(func, x0, vis.imag[i], p0)[0] )
	return npfmt(pf).flatten()


def _DoMultiprocess_Plot( iterable ) : 
	n1, n2 = iterable[0]
	array = iterable[1][:-1]
	strbl = iterable[1][-1]
	freq, func, color, outdir, which = iterable[2][:5]
	if (which == 'phaseadd') : bl = iterable[2][-1]
	for i in xrange(n2-n1) : 
		if (which == 'amp') : arraymax = np.array([0.,0])
		elif (which == 'phaseadd') : strl = '=(%.3f, %.3f)' % tuple(bl[i][:2])
		for j in xrange(len(array)) : 
			if (array[j] is None) : continue
			if (which == 'amp') : 
				nmf = array[j][i].size/5
				nmf = nmf if(nmf%2==1) else nmf+1
				arraymax[j] =spsn.medfilt(array[j][i],nmf).max()
			plt.plot(freq, array[j][i], color=color[j], ls='', marker='o', markersize=3, label=func[j])
		plt.xlim(freq.min(), freq.max())
		plt.xlabel(r'$\nu$ [MHz]', size=16)
		plt_axes('x', 'both', [10,2])
		if (which == 'amp') : 
			plt.ylabel('Amp [A.U.]', size=16)
			plt.ylim(0, 1.3*arraymax.min())
			plt.legend(loc=8)
			plt.title('Fitted amplitude of abs(visibility), baseline='+strbl[i], size=16)
			plt.savefig(outdir+'Amp_'+strbl[i]+'.png')
		elif (which == 'deff') : 
			plt.ylabel('Deff [m]', size=16)
			plt.ylim(4.5, 6.5)
			plt_axes('y', 'both', [0.1,0.05], '%.1f')
			plt.legend()
			plt.title('Fitted effective diameter, baseline='+strbl[i], size=16)
			plt.savefig(outdir+'Deff_'+strbl[i]+'.png')
		elif (which == 'phaseadd') : 
			plt.ylabel('Phaseadd [deg]', size=16)
			plt_axes('y', 'both', [30,5])
			plt.ylim(0, 360)
			plt.title('Fitted additional phase, baseline='+strbl[i]+strl, size=16)
			plt.legend()
			plt.savefig(outdir+'Phaseadd_'+strbl[i]+'.png')
		plt.close()


##################################################
##################################################
##################################################


class CaliPhase( object ) : 
 	dtype = 'class:'+sys._getframe().f_code.co_name


	def __init__( self, antarray=None, masking=None ) :
		''' antarray: instance of class:AntArray
		    masking:  instance of class:Masking '''
		if (antarray is not None) : 
			nhdf5, nhdf5transit = antarray.Hdf5.nhdf5, antarray.Hdf5.nhdf5transit[-1]
			if (nhdf5 != nhdf5transit) : Raise(Exception, 'antarray.Hdf5.nhdf5='+str(nhdf5)+' is NOT the file which contains the transit source (nhdf5='+str(nhdf5transit)+')')
			self.antarray = antarray
			self.antarray.SelectVisType()  # Ensure cross1
			self._Outdir()
		if (masking is not None) : self.masking = masking
		self.showprogress = False
		self.smoothdone = False


	def _Outdir( self ) : 
		outdir = sys.argv[0][:-3] + '_output/'
		if (outdir[:2]=='./') : outdir = outdir[2:]
		mkdir(outdir)
		indir = self.antarray.Hdf5.hdf5dir[:-1].split('/')[-1]+'/'
		outdir += indir 
		mkdir(outdir)
		self.outdir = outdir + self.dtype[6:]+'/'
		mkdir(self.outdir)


	def ShowProgress( self, showprogress=True ) : 
		self.showprogress = showprogress


	def RADec( self, RA=None, Dec=None ) : 
		''' Dec: Dec of the data/antenna pointing/calibration source
		RA: RA of the calibration source 
		angle in degree '''
		if (RA  is not None) : self.RA  = RA
		if (Dec is not None) : self.Dec = Dec


	def Fringe( self, nsigma=6, plotnv=None, plotfreq=None ) : 
		''' Take where is the fringe of the bright source
		nsigma: times of the sigma of the beam
		=> self.nhdf5, self.timerange, self.vis, self.timem '''
		inttime = self.antarray.Ant.inttime
		Nt = self.antarray.vis.shape[0]
		timerange1, timerange2, ncount = -10, Nt+10, 0
		while (timerange1<0 and timerange2>Nt) : 
			nsigma -= 0.1*ncount
			# Ideal Gaussian beam
			sigma = 300/self.antarray.Ant.freq.mean()/self.antarray.Ant.dishdiam/2.287
			sigma = Sph2Circ(sigma, self.Dec*np.pi/180)  # rad
			timerange = 24*3600/(2*np.pi)*nsigma*sigma/inttime # pix
			#--------------------------------------------------
			transittime = self.antarray.Hdf5.transittimelocal[-1] /inttime  # pixels
			timerange1, timerange2 = int(round(transittime-timerange/2)), int(round(transittime+timerange/2))  # pixels
		#--------------------------------------------------
		if (timerange2 > Nt) : 
			timerange1 = [timerange1, Nt]
			timerange2 = [0, timerange2-Nt]
			self.nhdf5 = (self.antarray.Hdf5.nhdf5, self.antarray.Hdf5.nhdf5+1)
			self.timerange = npfmt([timerange1, timerange2])
		elif (timerange1 < 0) : 
			timerange1 = [Nt+timerange1, Nt]
			timerange2 = [0, transittime2]
			self.nhdf5 = (self.antarray.Hdf5.nhdf5-1, self.antarray.Hdf5.nhdf5)
			self.timerange = npfmt([timerange1, timerange2])
		else : 
			self.nhdf5 = (self.antarray.Hdf5.nhdf5,)
			self.timerange = npfmt([[timerange1, timerange2]]) # pixels
		# self.timerange, self.nhdf5
		#--------------------------------------------------
		for i in xrange(len(self.nhdf5)) : 
			if (self.nhdf5[i] != self.antarray.Hdf5.nhdf5) : 
				antarray = AntArray(self.antarray.Hdf5.hdf5dir)
				antarray.WhichHdf5(self.nhdf5[i])
				antarray.MaskChannel(self.antarray.Blorder.maskchannel)
				antarray.SelectVisType(self.antarray.vistype)
				masking = Masking(antarray)
				try : masking.MaskNoiseSource(self.masking.noisesource.pixstart, self.masking.noisesource.pixlength, self.masking.noisesource.pixperiod)
				except : pass
				visno = antarray.vis[self.timerange[i,0]:self.timerange[i,1],:,self.antarray.visorder]
				maskno = masking.mask[self.timerange[i,0]:self.timerange[i,1]]
			else : 
				visyes = self.antarray.vis[self.timerange[i,0]:self.timerange[i,1],:,self.antarray.visorder]
				maskyes = self.masking.mask[self.timerange[i,0]:self.timerange[i,1]]
		#--------------------------------------------------
		if (len(self.nhdf5) == 1) : pass
		elif (self.nhdf5[0] == self.antarray.Hdf5.nhdf5) : 
			visyes = np.concatenate([visyes, visno], 0)
			maskyes = np.concatenate([maskyes, maskno], 0)
		else : 
			visyes = np.concatenate([visno, visyes], 0)
			maskyes = np.concatenate([maskno, maskyes], 0)
		self.vis = np.ma.MaskedArray(visyes, maskyes)
		self.timem = (np.arange(self.timerange[0,0], self.timerange[0,0]+self.vis.shape[0])+self.nhdf5[0]*Nt) *inttime/60. # min
		visyes = visno = maskyes = maskno = 0 #@
		if (self.antarray.visorder.size == 1) : 
			self.vis = self.vis[:,:,None]
		#--------------------------------------------------
		# Plot, longest East-West baseline, 750MHz
		if (plotnv is None) : 
			bl = self.antarray.Blorder.baseline[self.antarray.visorder]
			bl = abs(bl[:,0]) # longest East-West baseline
			nv = np.where(bl==bl.max())[0][0]
		else : nv = plotnv
		freq = self.antarray.Ant.freq
		if (plotfreq is None) : plotfreq = freq.mean()
		freqstr = str(int(round(plotfreq)))
		nvbl = np.arange(len(self.antarray.Blorder.blorder))[self.antarray.visorder][nv]
		bl = self.antarray.Blorder.Order2Bl(nvbl)
		strbl = str(bl[0])+'-'+str(bl[1])
		nf = abs(freq-plotfreq)
		nf = np.where(nf==nf.min())[0][0]
		vis = self.vis[:,nf,nv]
		vmax = abs(vis).max()*1.05
		plt.figure(figsize=(17,6))
		plt.subplot(1,2,1)
		plt.plot(self.timem, self.vis.real[:,nf,nv], 'b-')
		plt.xlabel('$t$ [min]', size=16)
		plt.xlim(self.timem.min(), self.timem.max())
		plt_axes('x', 'both', [10,1])
		plt.ylabel('[A.U.]', size=16)
		plt.ylim(-vmax, vmax)
		plt.title('Real part', size=16)
		plt.subplot(1,2,2)
		plt.plot(self.timem, self.vis.imag[:,nf,nv], 'r-')
		plt.xlabel('$t$ [min]', size=16)
		plt.xlim(self.timem.min(), self.timem.max())
		plt_axes('x', 'both', [10,1])
		plt.ylabel('[A.U.]', size=16)
		plt.ylim(-vmax, vmax)
		plt.title('Imaginary part', size=16)
		plt.suptitle('Fringe of baseline='+strbl+' @ '+freqstr+'MHz', fontsize=20)
		plt.savefig(self.outdir+'vis_'+strbl+'_'+freqstr+'MHz.png')
		plt.close()
		#--------------------------------------------------
		# Save files
	#	np.save(self.outdir+'vis.data_nf-'+str(nf)+'_nv-'+str(nv), vis.data)
	#	np.save(self.outdir+'vis.mask_nf-'+str(nf)+'_nv-'+str(nv), vis.mask, bool)
		np.save(self.outdir+'timem_minute', self.timem)



	def Smooth( self, timetimes=0, freqtimes=0 ) : 
#	def Smooth( self, timeper=0, timetimes=1, freqper=0, freqtimes=1 ) : 
		'''
		Fill/Reset the masked elements, and smooth.
		If just want to reset the masked elements, can set per=0
		'''
		vis = ResetMasked(self.vis, 0)  # return real array
		vis = Smooth(vis, 1, 3, freqtimes)
		vis = Smooth(vis, 0, 3, timetimes)
		self.vis = np.ma.MaskedArray(vis, self.vis.mask)
		self.smoothdone = True



	def FitBeam( self, Nprocess=None ) : 
		''' Fit the beam with abs(vis) 
		Nprocess: number of processes in multiprocessing '''
		Nprocess = NprocessCPU(Nprocess)[0]
		if (not self.smoothdone) : self.Smooth()
		# Initial guess
		inttime = self.antarray.Ant.inttime
		s0 = 300/self.antarray.Ant.freq.mean()/self.antarray.Ant.dishdiam/2.287
		t0 = self.vis.shape[0]/2*inttime  # second
		x0 = np.arange(self.vis.shape[0])*inttime
		#--------------------------------------------------
		beam = abs(self.vis.data)
		beam = ArrayAxis(beam, 0, -1, 'move')
		shape = beam.shape
		# Reshape to 2D
		beam = beam.reshape(np.prod(shape[:2]), shape[-1])
		#--------------------------------------------------
		bcast = (x0, t0, s0, self.Dec, self.showprogress)
		pool = PoolFor(0, len(beam), Nprocess)
		pf = pool.map_async(_DoMultiprocess_FitBeam, beam, bcast)
		#--------------------------------------------------
		pf = np.concatenate(pf).T
		pf = pf.reshape( (len(pf),) + shape[:2] )
		self.Ampb, self.Timeb, self.Sigmab = pf
		del pool, beam
		self.Ampb = abs(self.Ampb)
		self.Sigmab = abs(self.Sigmab)
		#--------------------------------------------------
		np.save(self.outdir+'Ampb.npy'  , self.Ampb  )
		np.save(self.outdir+'Timeb.npy' , self.Timeb )
		np.save(self.outdir+'Sigmab.npy', self.Sigmab)
		bl =self.antarray.Blorder.baseline[self.antarray.visorder]
		self.Phasens = 2*np.pi/300*self.antarray.Ant.freq[:,None]*(-bl[:,1][None,:])*np.sin((self.Dec-self.antarray.Ant.lonlat[1])*np.pi/180)
		np.save(self.outdir+'Phasens.npy', self.Phasens)
		# vis * exp(-1j*Phasens) * exp(-1j*Phaseadd)
		return [self.Ampb, self.Timeb, self.Sigmab, self.Phasens]



	def FitPhase( self, Nprocess=None ) : 
		Nprocess = NprocessCPU(Nprocess)[0]
	#	if (not self.smoothdone) : self.Smooth()
		if ('Timeb' not in self.__dict__.keys()) : self.FitBeam()
		bl =self.antarray.Blorder.baseline[self.antarray.visorder]
		freq = self.antarray.Ant.freq
		inttime = self.antarray.Ant.inttime
		x0 = np.arange(self.vis.shape[0])*inttime
	#	t0 = self.vis.shape[0]/2*inttime  # for smoothdone
		#--------------------------------------------------
		phase = np.angle(self.vis.data)
		phase = ArrayAxis(phase, 0, -1, 'move')
		shape = phase.shape
		# Reshape to 2D
		phase = phase.reshape(np.prod(shape[:2]), shape[-1])
		#--------------------------------------------------
		pool = PoolFor(0, len(phase), Nprocess)
	#	send = phase  # smoothdone
	#	bcast = (x0, t0, self.Dec, freq, bl, self.showprogress, shape[1], self.antarray.Ant.lonlat[1])  # smoothdone
		send = (phase, self.Timeb.flatten())
		bcast = (x0, self.Dec, freq, bl, self.showprogress, shape[1], self.antarray.Ant.lonlat[1])
		pf = pool.map_async(_DoMultiprocess_FitPhase, send, bcast)
		#--------------------------------------------------
		pf = np.concatenate(pf) %(2*np.pi)
		self.Phaseaddp = pf.reshape(shape[:2])  # (Nf,Nv)
		#--------------------------------------------------
		np.save(self.outdir+'Phaseaddp.npy', self.Phaseaddp)
		return self.Phaseaddp



	def FitVis(self, plotnv=None, Nprocess=None, plotfreq=None): 
		Nprocess = NprocessCPU(Nprocess)[0]
		if ('Timeb' not in self.__dict__.keys()) : self.FitBeam()
		bl =self.antarray.Blorder.baseline[self.antarray.visorder]
		freq = self.antarray.Ant.freq
		inttime = self.antarray.Ant.inttime
		s0 = 300/freq.mean()/self.antarray.Ant.dishdiam/2.287
		A0 = abs(self.vis.data).max(0).flatten()
		x0 = np.arange(self.vis.shape[0])*inttime  # second
		#--------------------------------------------------
		# Move time axis to -1
		vis = ArrayAxis(self.vis.data.imag, 0, -1, 'move')
	#	vis = ArrayAxis(self.vis.data.real, 0, -1, 'move')
		shape = vis.shape  # (512,528,3600)
		# Reshape to 2D
		vis = vis.reshape(np.prod(shape[:2]), shape[-1])
		#--------------------------------------------------
		pool = PoolFor(0, len(vis), Nprocess)
		send = (vis, A0, self.Timeb.flatten())
		bcast = (x0, s0, self.Dec, freq, bl, self.showprogress, shape[1], self.antarray.Ant.lonlat[1])
		pf = pool.map_async(_DoMultiprocess_FitVis, send, bcast)
		#--------------------------------------------------
		pf = np.concatenate(pf).T  # shape=(4, Nf*Nv)
		del pool, bcast, vis, A0
		pf = pf.reshape( (len(pf),) + shape[:2] )
		self.Ampv, self.Sigmav, self.Phaseaddv = pf
		self.Sigmav = abs(self.Sigmav)
		self.Phaseaddv[self.Ampv<0] += np.pi
		self.Phaseaddv %= (2*np.pi)
		self.Ampv = abs(self.Ampv)
		#--------------------------------------------------
		np.save(self.outdir+'Ampv.npy'     , self.Ampv     )
		np.save(self.outdir+'Sigmav.npy'   , self.Sigmav   )
		np.save(self.outdir+'Phaseaddv.npy', self.Phaseaddv)
		#--------------------------------------------------
		# Plot
		if (plotnv is None) : 
			bl = abs(bl[:,0]) # longest East-West baseline
			nv = np.where(bl==bl.max())[0][0]
		else : nv = plotnv
		nvbl = np.arange(len(self.antarray.Blorder.blorder))[self.antarray.visorder][nv]
		bl = self.antarray.Blorder.Order2Bl(nvbl)
		strbl = str(bl[0])+'-'+str(bl[1])
		bl =self.antarray.Blorder.baseline[self.antarray.visorder]
		if (plotfreq is None) : plotfreq = freq.mean()
		nf = abs(freq-plotfreq)
		freqstr = str(int(round(plotfreq)))
		nf = np.where(nf==nf.min())[0][0]
		vis = self.vis[:,nf,nv]  # real+imag+masked
		vmax = abs(vis).max()*1.05
		#--------------------
		def func(x, p) : 
			beta = Circ2Sph((x-p[1])/3600.*15*np.pi/180, self.Dec*np.pi/180)
			A = p[0] * np.exp(-beta**2 /2 /p[2]**2)
			angle = 2*np.pi/(300./freq[nf]) * ( bl[nv,0]*np.sin(beta) - bl[nv,1]*np.cos(beta)*np.sin((self.Dec-self.antarray.Ant.lonlat[1])*np.pi/180) ) + p[3]
			return A*np.exp(1j*angle)  # complex
		#--------------------
		visf = func(x0, [self.Ampv[nf,nv], self.Timeb[nf,nv], self.Sigmav[nf,nv], self.Phaseaddv[nf,nv]])
		plt.figure(figsize=(17,6))
		plt.subplot(1,2,1)
		plt.plot(self.timem, vis.real, 'b-', label='data')
		plt.plot(self.timem, visf.real, 'r-', label='Fitting')
		plt.legend(fontsize=16)
		plt.xlabel('$t$ [min]', size=16)
		plt.xlim(self.timem.min(), self.timem.max())
		plt_axes('x', 'both', [10,1])
		plt.ylabel('[A.U.]', size=16)
		plt.ylim(-vmax, vmax)
		plt.title('Real part', size=16)
		plt.subplot(1,2,2)
		plt.plot(self.timem, vis.imag, 'b-', label='data')
		plt.plot(self.timem, visf.imag, 'r-', label='Fitting')
		plt.legend(fontsize=16)
		plt.xlabel('$t$ [min]', size=16)
		plt.xlim(self.timem.min(), self.timem.max())
		plt_axes('x', 'both', [10,1])
		plt.ylabel('[A.U.]', size=16)
		plt.ylim(-vmax, vmax)
		plt.title('Imaginary part', size=16)
		plt.suptitle('Fitted fringe of baseline='+strbl+' @ '+freqstr+'MHz by FitVis()', fontsize=20)
		plt.savefig(self.outdir+'vis_FitVis_'+strbl+'_'+freqstr+'MHz.png')
		plt.close()
		return [self.Ampv, self.Sigmav, self.Phaseaddv]



	def FitExp( self, plotnv=None, Nprocess=None, plotfreq=None):
		Nprocess = NprocessCPU(Nprocess)[0]
		if ('Timeb' not in self.__dict__.keys()) : self.FitBeam()
		bl =self.antarray.Blorder.baseline[self.antarray.visorder]
		freq = self.antarray.Ant.freq
		inttime = self.antarray.Ant.inttime
		x0 = np.arange(self.vis.shape[0])*inttime  # second
		#--------------------------------------------------
		# Move time axis to -1
		vis = ArrayAxis(self.vis.data, 0, -1, 'move')
		shape = vis.shape  # (512,528,3600)
		# Reshape to 2D
		vis = vis.reshape(np.prod(shape[:2]), shape[-1])
		#--------------------------------------------------
		pool = PoolFor(0, len(vis), Nprocess)
		send = (vis, self.Timeb.flatten())
		bcast = (x0, self.Dec, freq, bl, self.showprogress, shape[1], self.antarray.Ant.lonlat[1])
		pf = pool.map_async(_DoMultiprocess_FitExp, send, bcast)
		#--------------------------------------------------
	#	self.Phaseadde = np.concatenate(pf) %(2*np.pi)
	#	self.Phaseadde = self.Phaseadde.reshape(shape[:2])
		pf = np.concatenate(pf).reshape(shape[:2])
		self.Phaseadde = pf %(2*np.pi)
		del pool, bcast, vis
		#--------------------------------------------------
		np.save(self.outdir+'Phaseadde.npy', self.Phaseadde)
		#--------------------------------------------------
		# Plot
		if (plotnv is None) : 
			bl = abs(bl[:,0]) # longest East-West baseline
			nv = np.where(bl==bl.max())[0][0]
		else : nv = plotnv
		nvbl = np.arange(len(self.antarray.Blorder.blorder))[self.antarray.visorder][nv]
		bl = self.antarray.Blorder.Order2Bl(nvbl)
		strbl = str(bl[0])+'-'+str(bl[1])

		bl =self.antarray.Blorder.baseline[self.antarray.visorder]
		if (plotfreq is None) : plotfreq = freq.mean()
		freqstr = str(int(round(plotfreq)))
		nf = abs(freq-plotfreq)
		nf = np.where(nf==nf.min())[0][0]
		vis = self.vis[:,nf,nv]  # real+imag+masked
		vmax = abs(vis).max()*1.05
		A = abs(vis.data)
		#--------------------
		def func(x, p) : 
			beta = Circ2Sph((x-self.Timeb[nf,nv])/3600.*15*np.pi/180, self.Dec*np.pi/180)
			angle = 2*np.pi/(300./freq[nf]) * ( bl[nv,0]*np.sin(beta) - bl[nv,1]*np.cos(beta)*np.sin((self.Dec-self.antarray.Ant.lonlat[1])*np.pi/180) ) + p[0]
			return A*np.exp(1j*angle)  # complex
		#--------------------
		visf = func(x0, [self.Phaseadde[nf,nv]])
		plt.figure(figsize=(17,6))
		plt.subplot(1,2,1)
		plt.plot(self.timem, vis.real, 'b-', label='data')
		plt.plot(self.timem, visf.real, 'r-', label='Fitting')
		plt.legend(fontsize=16)
		plt.xlabel('$t$ [min]', size=16)
		plt.xlim(self.timem.min(), self.timem.max())
		plt_axes('x', 'both', [10,1])
		plt.ylabel('[A.U.]', size=16)
		plt.ylim(-vmax, vmax)
		plt.title('Real part', size=16)
		plt.subplot(1,2,2)
		plt.plot(self.timem, vis.imag, 'b-', label='data')
		plt.plot(self.timem, visf.imag, 'r-', label='Fitting')
		plt.legend(fontsize=16)
		plt.xlabel('$t$ [min]', size=16)
		plt.xlim(self.timem.min(), self.timem.max())
		plt_axes('x', 'both', [10,1])
		plt.ylabel('[A.U.]', size=16)
		plt.ylim(-vmax, vmax)
		plt.title('Imaginary part', size=16)
		plt.suptitle('Fitted fringe of baseline='+strbl+' @ '+freqstr+'MHz by FitExp()', fontsize=20)
		plt.savefig(self.outdir+'vis_FitExp_'+strbl+'_'+freqstr+'MHz.png')
		plt.close()
		return self.Phaseadde



	def Plot( self, Nprocess=None ) : 
		outdir = self.outdir + 'Figure/'
		mkdir(outdir)
		qt = ['amp','deff','phaseadd']
		#--------------------------------------------------
		nbl = np.arange(len(self.antarray.visorder))
		bl =self.antarray.Blorder.blorder[self.antarray.visorder]
		strbl = []
		for i in xrange(len(nbl)) : 
			strbl.append(str(bl[i,0])+'-'+str(bl[i,1]))
		#--------------------------------------------------
		self._PlotAmp(outdir, strbl, Nprocess)
		Raise()
		self._PlotDeff(outdir, strbl, Nprocess)
		self._PlotPhaseadd(outdir, strbl, Nprocess)


	def _PlotAmp( self, outdir, strbl, Nprocess ) : 
		'''
		whichbl: 
			1D integer ndarray/list/tuple: visorder of the baseline
			2D integer ndarray(shape=(N,2)) / list[(1,3),(1,7),...] / tuple((1,3),(2,8),...): channel pair
		'''
		freq = self.antarray.Ant.freq
		func, amp, color = ['FitBeam', 'FitVis'], [], ['b','r']
		try : amp.append(self.Ampb.T)  # (Nv,Nf)
		except : amp.append(None)
		try : amp.append(self.Ampv.T)
		except : amp.append(None)
		if (amp[0] is None and amp[1] is None): 
			Raise(Warning, "CaliPhase.Plot('Amp'): NOT exist self.Ampb, self.Ampv")
			return
		send = tuple(amp) + (strbl,)
		bcast = (freq, func, color, outdir, 'amp')
		print 'enenen amp'
		Raise()
		pool = PoolFor(0, len(strbl), Nprocess)
		pool.map_async(_DoMultiprocess_Plot, send, bcast)


	def _PlotDeff( self, outdir, strbl, Nprocess ) : 
		freq = self.antarray.Ant.freq
		func, Deff, color = ['FitBeam', 'FitVis'], [], ['b','r']
		try : Deff.append(300/(self.Sigmab.T*2.287*freq[None,:]))
		except : Deff.append(None)
		try : Deff.append(300/(self.Sigmav.T*2.287*freq[None,:]))
		except : Deff.append(None)
		if (Deff[0] is None and Deff[1] is None) : 
			Raise(Warning, "CaliPhase.Plot('Deff'): NOT exist self.Sigmab, self.Sigmav")
			return
		send = tuple(Deff) + (strbl,)
		bcast = (freq, func, color, outdir, 'deff')
		pool = PoolFor(0, len(strbl), Nprocess)
		pool.map_async(_DoMultiprocess_Plot, send, bcast)


	def _PlotPhaseadd( self, outdir, strbl, Nprocess ) : 
		freq = self.antarray.Ant.freq
	#	func,phaseadd,color = ['FitPhase','FitExp'], [], ['b','r']
		func,phaseadd,color = ['FitPhase','FitVis','FitExp'], [], ['b','r','c']
		try : phaseadd.append(self.Phaseaddp.T*180/np.pi) # deg
		except : phaseadd.append(None)
		try : phaseadd.append(self.Phaseaddv.T*180/np.pi)
		except : phaseadd.append(None)
		try : phaseadd.append(self.Phaseadde.T*180/np.pi)
		except : phaseadd.append(None)
		if (phaseadd[0] is None and phaseadd[1] is None and phaseadd[2] is None) : 
			Raise(Warning, "CaliPhase.Plot('Phaseadd'): NOT exist self.Phaseaddp, self.Phaseaddv, self.Phaseadde")
			return
		bl =self.antarray.Blorder.baseline[self.antarray.visorder]
		send = tuple(phaseadd) + (strbl,)
		bcast = (freq, func, color, outdir, 'phaseadd', bl)
		pool = PoolFor(0, len(strbl), Nprocess)
		pool.map_async(_DoMultiprocess_Plot, send, bcast)

