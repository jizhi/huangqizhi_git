import os
import time
from AntArray import *
from Masking import *
import jizhipy as jp
from jizhipy.Plot import *

#import scipy.signal as spsn
#from Sph2Circ import *
#from Plot import *
#from PoolFor import *
#from ResetMasked import *
#from Smooth import *
#from FuncFit import *

##################################################



def _DoMultiprocess_FitBeam( iterable ) : 
	n1, n2  = iterable[0]
	beam = iterable[1]
	timem, s0, Dec, verbose = iterable[2]
	if (verbose) : progressbar = jp.ProgressBar('    pid='+str(os.getpid())+':', len(beam), False)
	#--------------------------------------------------
	def func(x, p) :  # x is time, not angle
		beta = jp.Circ2Sph((x-p[1])/60.*15*np.pi/180, Dec*np.pi/180)
		y = p[0] * np.exp(-beta**2 /2 /p[2]**2)
		return y
	#--------------------------------------------------
	pf = []
	for i in xrange(n2-n1) : 
		if (verbose) : progressbar.Progress()
		p0 = [beam[i].max(), timem.mean(), s0]
		pf.append( jp.FuncFit(func, timem, beam[i], p0, warning=False)[0] )
	return jp.npfmt(pf)



def _DoMultiprocess_FitPhase( iterable ) : 
	n1, n2 = iterable[0]
	phase, t0list = iterable[1]
	timem, Dec, freq, bl, verbose, Nf, lat = iterable[2]
	if (verbose) : progressbar = jp.ProgressBar('    pid='+str(os.getpid())+':', len(phase), False)
	pf, plf = [], []
	#--------------------------------------------------
	for i in xrange(n2-n1) : 
		t0 = t0list[i] #@
		if (verbose) : progressbar.Progress()
		nv, nf = (i+n1)/Nf, (i+n1)%Nf
		#--------------------------------------------------
		def funcNO(x, p) : 
			beta =jp.Circ2Sph((x-t0)/60.*15*np.pi/180, Dec*np.pi/180)
			angle = 2*np.pi/(300./freq[nf]) *( bl[nv,0]*np.sin(beta) - bl[nv,1]*np.cos(beta)*np.sin((Dec-lat)*np.pi/180) ) + p[0]
			return np.sin(angle)
		pf.append( list(jp.FuncFit(funcNO, timem, np.sin(phase[i]), [1])[0]) )
		#--------------------------------------------------
#		def funcYES(x, p) : 
#			beta =jp.Circ2Sph((x-t0)/60.*15*np.pi/180, Dec*np.pi/180)
#			angle = 2*np.pi/(300./freq[nf]) *( p[0]*np.sin(beta) - bl[nv,1]*np.cos(beta)*np.sin((Dec-lat)*np.pi/180) ) + pf[-1][0]
#			return np.sin(angle)
#		pf[-1] += list(jp.FuncFit(funcYES, timem, np.sin(phase[i]), bl[nv,0])[0])
		#--------------------------------------------------
		def funcBOTH(x, p) : 
			beta =jp.Circ2Sph((x-t0)/60.*15*np.pi/180, Dec*np.pi/180)
			angle = 2*np.pi/(300./freq[nf]) *( p[1]*np.sin(beta) - bl[nv,1]*np.cos(beta)*np.sin((Dec-lat)*np.pi/180) ) + p[0]
			return np.sin(angle)
		plf.append( list(jp.FuncFit(funcBOTH, timem, np.sin(phase[i]), [pf[-1][0], bl[nv,0]])[0]) )
#	return jp.npfmt(pf)
	return jp.npfmt(plf)



def _DoMultiprocess_FitVis( iterable ) : 
	n1, n2 = iterable[0]
	vis, A0, t0list = iterable[1]
	x0, s0, Dec, freq, bl, showprogress, Nv, lat = iterable[2]
	if (showprogress) : progressbar = jp.ProgressBar('CaliPhase.FitPhase(), pid='+str(os.getpid())+':', len(phase))
	pf = []
	#--------------------------------------------------
	for i in xrange(n2-n1) : 
		if (showprogress) : progressbar.Progress()
		nf, nv = (i+n1)/Nv, (i+n1)%Nv
		def func(x, p) : 
			beta = jp.Circ2Sph((x-t0list[i])/3600.*15*np.pi/180, Dec*np.pi/180)
			A = p[0] * np.exp(-beta**2 /2 /p[1]**2)
			angle = 2*np.pi/(300./freq[nf]) *( bl[nv,0]*np.sin(beta) - bl[nv,1]*np.cos(beta)*np.sin((Dec-lat)*np.pi/180) ) + p[2]
			return A*np.sin(angle)
		#	return A*np.cos(angle)
		#--------------------------------------------------
		p0 = [A0[i], s0, 1]
		pf.append( jp.FuncFit(func, x0, vis[i], p0)[0] )
	return jp.npfmt(pf)  # shape=(n2-n1, 4)


def _DoMultiprocess_FitExp( iterable ) : 
	n1, n2 = iterable[0]
	vis, t0list = iterable[1]
	x0, Dec, freq, bl, showprogress, Nv, lat = iterable[2]
	if (showprogress) : progressbar = jp.ProgressBar('CaliPhase.FitExp(), pid='+str(os.getpid())+':', len(phase))
	pf = []
	#--------------------------------------------------
	for i in xrange(n2-n1) : 
		if (showprogress) : progressbar.Progress()
		nf, nv = (i+n1)/Nv, (i+n1)%Nv
		beta = jp.Circ2Sph((x0-t0list[i])/3600.*15*np.pi/180, Dec*np.pi/180)
		A = abs(vis[i])
		def func(x, p) : 
			angle = 2*np.pi/(300./freq[nf]) *( bl[nv,0]*np.sin(beta) - bl[nv,1]*np.cos(beta)*np.sin((Dec-lat)*np.pi/180) ) + p[0]
			return A*np.sin(angle)
		#	return A*np.cos(angle)
		#--------------------------------------------------
		p0 = np.pi
		pf.append( jp.FuncFit(func, x0, vis.imag[i], p0)[0] )
	return jp.npfmt(pf).flatten()


def _DoMultiprocess_Plot( iterable ) : 
	n1, n2 = iterable[0]
	array = iterable[1][:-1]
	strbl = iterable[1][-1]
	freq, func, color, outdir, which = iterable[2][:5]
	which = which.lower()
	if (which == 'deff') : 
		dyDeff = iterable[2][5]
		fwhm2sigmafactor = iterable[2][6]
		dishdiam = iterable[2][7]
	elif (which in ['lew', 'phaseadd', 'phase']) : 
		bl = iterable[2][5]
		if (which == 'lew') : dyLew = iterable[2][6]
	#--------------------------------------------------
	for i in xrange(n2-n1) : # baseline
		if (which == 'amp') : arraymax = np.array([0.,0])+1e30
		elif (which in ['lew', 'phaseadd', 'phase') : strl = '=(%.3f, %.3f)' % tuple(bl[i][:2])
		for j in xrange(len(array)) : # FitBeam, FitVis, ...
			if (array[j] is None) : continue
			nmf = array[j][i].size/5
			nmf = nmf if(nmf%2==1) else nmf+1
			if (which == 'amp') : 
				arraymax[j] = spsn.medfilt(array[j][i], nmf).max()
			plt.plot(freq, array[j][i], color=color[j], ls='', marker='o', markersize=3, label=func[j])
		plt.xlim(int(round(freq.min())), int(round(freq.max())))
		plt.xlabel(r'$\nu$ [MHz]', size=16)
		plt_axes('x', 'both', [25,5])
		#--------------------------------------------------
		if (which == 'amp') : 
			plt.ylabel('Amp [A.U.]', size=16)
			plt.ylim(0, 1.4*arraymax.min())
		#	plt_axisformat('y')
			plt.legend()
			plt.title('Fitted amplitude of abs(visibility), baseline='+strbl[i], size=16)
			plt.savefig(outdir+'Amp_'+strbl[i]+'.png')
		#--------------------------------------------------
		elif (which == 'deff') : 
			plt.ylabel(r'Deff$=\frac{\lambda}{%.3f \cdot \sigma}$ [m]' % fwhm2sigmafactor, size=16)
			if (dyDeff is not None) : 
				amean = spsn.medfilt(array[j][i], nmf).mean()
			#	plt.ylim(amean-dyDeff, dishdiam*1.05)
				plt.ylim(amean-dyDeff, amean+dyDeff)
			plt_axes('y', 'both', [0.1, 0.05], '%.1f')
			plt.legend()
		#	plt.title(r'Fitted effective diameter Deff$=\frac{\lambda}{'+('%.3f' % fwhm2sigmafactor)+r' \cdot \sigma}$, baseline='+strbl[i], size=16)
			plt.title(r'Fitted effective diameter ('+str(dishdiam)+'m), baseline='+strbl[i], size=16)
			plt.savefig(outdir+'Deff_'+strbl[i]+'.png')
		#--------------------------------------------------
		elif (which == 'lew') : 
			plt.ylabel(r'$L_{ew}$ [m]', size=16)
			if (dyLew is not None) : 
				amean = spsn.medfilt(array[j][i], nmf).mean()
			#	sign = np.sign(amean[i][0])
			#	if (sign == 0) : sign = 1
			#	a1, a2 = abs(amean)-dyLew, abs(amean)+dyLew
			#	b1, b2 = abs(bl[i][0]-dyLew), abs(bl[i][0]+dyLew)
				plt.ylim(amean-dyLew, amean+dyLew)
			plt_axes('y', 'both', [0.05, 0.01], '%.1f')
			plt.legend()
			plt.title(r'Fitted Lew, baseline='+strbl[i]+strl,size=16)
			plt.savefig(outdir+'Lew_'+strbl[i]+'.png')
		#--------------------------------------------------
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


	def __init__( self, antarray=None, masking=None, fwhm2sigmafactor=2.287, Nprocess=None, verbose=True ) : 
		'''
		antarray: instance of class:AntArray
		masking:  instance of class:Masking
		caligain: instance of class:CaliGain

		Nprocess=100 is fast: 17 second
		'''
		self.verbose = verbose
		try : self.fwhm2sigmafactor = float(fwhm2sigmafactor)
		except : self.fwhm2sigmafactor = 2.287
		self.Nprocess = jp.NprocessCPU(Nprocess, verbose)[0]
		if (antarray is not None) : 
			if (antarray.vistype not in ['cross1','cross2']) : 
				jp.Raise(Exception, 'antarray.vistype="'+antarray.vistype+'" not in ["cross1", "cross2"]')
			nhdf5, nhdf5transit = antarray.Hdf5.nhdf5, antarray.Hdf5.nhdf5transit[-1]
			if (nhdf5 != nhdf5transit) : jp.Raise(Exception, 'antarray.Hdf5.nhdf5='+str(nhdf5)+' is NOT the file which contains the transit source (nhdf5='+str(nhdf5transit)+')')
			self.antarray = antarray  # share address, not copy
		#	self.antarray.SelectVisType()  # Ensure cross1
		if (masking is not None) : self.masking = masking
	#	if (caligain is not None) : self.caligain = caligain
		self.outdir = jp.Outdir((None,'file'), (0,'file'))



	def RADec( self, RA=None, Dec=None ) : 
		'''
		Dec: Dec of the data/antenna pointing/calibration source
		RA : RA  of the calibration source 
		angle in degree
		'''
		if (RA  is not None) : self.RA  = RA
		if (Dec is not None) : self.Dec = Dec



	def Fringe( self, nsigma=4, plotnv=None, plotfreq=None ) : 
		'''
		Take where is the fringe of the bright source
		nsigma: times of the sigma of the beam
		Generate: 
			self.nhdf5: nhdf5 of the fringe
			self.timerange: pixels range of the fringe
			self.vis: MaskedArray of fringe data, has been reset by maskvalue
			self.timem: in minute, self.timem.size==self.vis.shape[0]
			self.shift: a constant shift


		If there is an extra constant input/interference in front of A/D, this constant will become delta function (real value) after FFT, that will shift the real part of visibilities by a constant, but doesn't affect the imaginary part.

		For auto-correlation, we can't remove this shifting because we don't know how much it is, however, for cross-correlation, the base line/average? should be zero because of exp{ix}, we can remove the shifting in the cross-correlation.
		'''
		if (self.verbose) : print 'CaliPhase.Fringe: start @', jp.Time(1)
		self.shift = []
		inttime = self.antarray.Ant.inttime
		freq = self.antarray.Ant.freq
		Nt = self.antarray.vis.shape[0]  # each hdf5 must have the same number of points
		timerange1, timerange2, ncount = -10, Nt+10, 0
		# Here we assume: fringe must be in 1 or 2 hdf5 files, not more than 2
		#--------------------------------------------------
		while (timerange1<0 and timerange2>Nt) : 
			nsigma -= 0.1*ncount
			# Ideal Gaussian beam
			sigma = 300/freq.mean()/self.antarray.Ant.dishdiam/self.fwhm2sigmafactor
			sigma = jp.Sph2Circ(sigma, self.Dec*np.pi/180)  # rad
			timerange = 24*3600/(2*np.pi)*nsigma*sigma/inttime # total pixels of the fringe in nsigma
			# transittime is calculated by LAST, starts from 0
			transittime = self.antarray.Hdf5.transittimelocal[-1] /inttime  # pixels
			timerange1, timerange2 = int(round(transittime-timerange/2)), int(round(transittime+timerange/2))  # pixels
		#--------------------------------------------------
		if (timerange2 > Nt) : 
			timerange1 = [timerange1, Nt]
			timerange2 = [0, timerange2-Nt]
			self.nhdf5 = (self.antarray.Hdf5.nhdf5, self.antarray.Hdf5.nhdf5+1)
			self.timerange = jp.npfmt([timerange1, timerange2])
		elif (timerange1 < 0) : 
			timerange1 = [Nt+timerange1, Nt]
			timerange2 = [0, transittime2]
			self.nhdf5 = (self.antarray.Hdf5.nhdf5-1, self.antarray.Hdf5.nhdf5)
			self.timerange = jp.npfmt([timerange1, timerange2])
		else : 
			self.nhdf5 = (self.antarray.Hdf5.nhdf5,)
			# self.timerange is range of pixels
			self.timerange = jp.npfmt([[timerange1, timerange2]])
		# Above generate: self.timerange, self.nhdf5
		#--------------------------------------------------
		for i in xrange(len(self.nhdf5)) : 
			if (self.nhdf5[i] != self.antarray.Hdf5.nhdf5) : 
				# Read self.nhdf5[i]
				antarray = AntArray(self.antarray.Hdf5.hdf5dir)
				antarray.WhichHdf5(self.nhdf5[i])
				antarray.MaskChannel(self.antarray.Blorder.maskchannel)
				antarray.SelectVisType(self.antarray.vistype)
				# mask it
				masking = Masking(antarray)
				try : masking.MaskNoiseSource(self.masking.noisesource.pixstart, self.masking.noisesource.pixlength, self.masking.noisesource.pixperiod)
				except : pass
				try : 
					n = len(self.Params.MaskLoopParams.__dict__.keys())
					for j in xrange(n) : 
						params = masking.Params.MaskLoopParams.__dict__['p'+str(j+1)]
						masking.MaskLoop(params.axis, params.per, params.times, params.nsigma, params.nloop, params.threshold)
				except : pass
				visno = antarray.vis[:,:,self.antarray.visorder]
				if (len(visno.shape) == 2) : visno = visno[:,:,None]
				try : visno[masking.mask] = maskvalue
				except : pass
				self.shift.append(visno.mean(0))
				visno = visno[self.timerange[i,0]:self.timerange[i,1]]
				maskno = masking.mask[self.timerange[i,0]:self.timerange[i,1]]
			else : 
				visyes = self.antarray.vis[:,:,self.antarray.visorder]
				if (len(visyes.shape)==2) : visyes = visyes[:,:,None]
				try : visyes[self.masking.mask] = self.masking.maskvalue
				except : pass
				self.shift.append(visyes.mean(0))
				visyes=visyes[self.timerange[i,0]:self.timerange[i,1]]
				maskyes = self.masking.mask[self.timerange[i,0]:self.timerange[i,1]]
		#--------------------------------------------------
		if (len(self.shift) == 2) : self.shift = ((self.shift[0]+self.shift[1])/2.)[None,:,:]
		else : self.shift = self.shift[0][None,:,:]
		#--------------------------------------------------
		if (len(self.nhdf5) == 1) : pass
		elif (self.nhdf5[0] == self.antarray.Hdf5.nhdf5) : 
			visyes = np.concatenate([visyes, visno], 0)
			maskyes = np.concatenate([maskyes, maskno], 0)
		else : 
			visyes = np.concatenate([visno, visyes], 0)
			maskyes = np.concatenate([maskno, maskyes], 0)
		self.vis = np.ma.MaskedArray(visyes, maskyes)  #@#@
		self.vis -= self.shift
		visyes = maskyes = visno = maskno = 0 #@
		#--------------------------------------------------
		self.timem = (np.arange(self.timerange[0,0], self.timerange[0,0]+self.vis.shape[0])+self.nhdf5[0]*Nt) *inttime/60. #@#@ min
		visyes = visno = maskyes = maskno = 0 #@
		#--------------------------------------------------
		#--------------------- Plot -----------------------
		np.save(self.outdir+'timem_minute', self.timem)
		self._PlotFringe(self.vis, ['Real part', 'Imaginary part'], ['b', 'r'])
		if (self.verbose) : print 'CaliPhase.Fringe:  end  @', jp.Time(1)+'\n'



	def Smooth( self, timetimes=100 ) : 
		''' Smooth(self.vis, 0, 3, timetimes) '''
		if (self.verbose) : print 'CaliPhase.Smooth: start @', jp.Time(1)
		array = [self.vis]
		mask = self.vis.mask.copy()
		self.vis = jp.Smooth(self.vis.data, 0, 3, timetimes, Nprocess=self.Nprocess)
		self.vis = np.ma.MaskedArray(self.vis, mask)
		array.append(self.vis)
		self._PlotFringe(array, [['Real part','Imaginary part'],['Smoothed','Smoothed']], [['b','r'],['m','g']], [1,3])
		if (self.verbose) : print 'CaliPhase.Smooth:  end  @', jp.Time(1)+'\n'



	def FitBeam( self ) : 
		''' Fit the beam with abs(vis) '''
		if (self.verbose) : print 'CaliPhase.FitBeam: start @', jp.Time(1)
		# Initial guess
		inttime = self.antarray.Ant.inttime
		s0 = 300/self.antarray.Ant.freq.mean()/self.antarray.Ant.dishdiam/2.287
		#--------------------------------------------------
		beam = abs(self.vis.data)
	#	beam = ArrayAxis(beam, 0, -1, 'move')
		beam = beam.T  # Faster
		shape = beam.shape
		# Reshape to 2D
		beam = beam.reshape(np.prod(shape[:2]), shape[-1])
		#--------------------------------------------------
		bcast = (self.timem, s0, self.Dec, self.verbose)
		pool = jp.PoolFor(0, len(beam), self.Nprocess)
		pf = pool.map_async(_DoMultiprocess_FitBeam, beam, bcast)
		print
		#--------------------------------------------------
		pf = np.concatenate(pf)  # (nv*nf, 3)
		pf = pf.reshape(shape[:2]+(pf.shape[-1],)).T  # (3,nf,nv)
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
		if (self.verbose) : print 'CaliPhase.FitBeam:  end  @', jp.Time(1)+'\n'



	def FitPhase( self ) : 
		'''
		(1) Because there is dPhase, fringe is not sensitive enough to Lew, we need to first fit the initial guess of dPhase.
		(2) Also because fringe is not sensitive enough to Lew, different range of fringe (nsigma=1, 2, 3, ...) have very different Lew.
		'''
		if ('Timeb' not in self.__dict__.keys()) : self.FitBeam()
		if (self.verbose) : print 'CaliPhase.FitPhase: start @', jp.Time(1)
		bl =self.antarray.Blorder.baseline[self.antarray.visorder]
		freq = self.antarray.Ant.freq
		inttime = self.antarray.Ant.inttime
		#--------------------------------------------------
		phase = np.angle(self.vis.data)
	#	phase = ArrayAxis(phase, 0, -1, 'move')
		phase = phase.T
		shape = phase.shape
		# Reshape to 2D
		phase = phase.reshape(np.prod(shape[:2]), shape[-1])
		#--------------------------------------------------
		pool = jp.PoolFor(0, len(phase), self.Nprocess)
		send = (phase, self.Timeb.T.flatten())
		bcast = (self.timem, self.Dec, freq, bl, self.verbose, shape[1], self.antarray.Ant.lonlat[1])
		pf = pool.map_async(_DoMultiprocess_FitPhase, send, bcast)
		print
		#--------------------------------------------------
		pf = np.concatenate(pf)  # (nv*nf, 2)
		pf = pf.reshape(shape[:2]+(pf.shape[-1],)).T  # (2,nf,nv)
		self.Phaseaddp, self.Lewp = pf
		self.Phaseaddp %= (2*np.pi)
		#--------------------------------------------------
		np.save(self.outdir+'Lewp.npy', self.Lewp)
		np.save(self.outdir+'Phaseaddp.npy', self.Phaseaddp)
		if (self.verbose) : print 'CaliPhase.FitPhase:  end  @', jp.Time(1)+'\n'



	def FitVis(self, plotnv=None, Nprocess=None, plotfreq=None): 
		if ('Timeb' not in self.__dict__.keys()) : self.FitBeam()
		bl =self.antarray.Blorder.baseline[self.antarray.visorder]
		freq = self.antarray.Ant.freq
		inttime = self.antarray.Ant.inttime
		s0 = 300/freq.mean()/self.antarray.Ant.dishdiam/2.287
		A0 = abs(self.vis.data).max(0).flatten()
		x0 = np.arange(self.vis.shape[0])*inttime  # second
		#--------------------------------------------------
		# Move time axis to -1
	#	vis = ArrayAxis(self.vis.data.imag, 0, -1, 'move')
		vis = vis.data.imag.T
		shape = vis.shape  # (512,528,3600)
		# Reshape to 2D
		vis = vis.reshape(np.prod(shape[:2]), shape[-1])
		#--------------------------------------------------
		pool = jp.PoolFor(0, len(vis), Nprocess)
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
			beta = jp.Circ2Sph((x-p[1])/3600.*15*np.pi/180, self.Dec*np.pi/180)
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
		pool = jp.PoolFor(0, len(vis), Nprocess)
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
			beta = jp.Circ2Sph((x-self.Timeb[nf,nv])/3600.*15*np.pi/180, self.Dec*np.pi/180)
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



	def Plot( self, dyDeff=None, dyLew=None, Nprocess=None ) : 
		if (self.verbose) : print 'CaliPhase.Plot: start @', jp.Time(1)
		outdir = jp.Mkdir(self.outdir + 'Plot/')
		#--------------------------------------------------
		nbl = np.arange(len(self.antarray.visorder))
		bl =self.antarray.Blorder.blorder[self.antarray.visorder]
		strbl = []
		for i in xrange(len(nbl)) : 
			strbl.append(str(bl[i,0])+'-'+str(bl[i,1]))
		#--------------------------------------------------
		self._PlotAmp( outdir, strbl, Nprocess)
		self._PlotDeff(outdir, strbl, Nprocess, dyDeff)
		self._PlotLew( outdir, strbl, Nprocess, dyLew)
		self._PlotPhaseadd(outdir, strbl, Nprocess)
		if (self.verbose) : print 'CaliPhase.Plot:  end  @', jp.Time(1)+'\n'



	def _PlotAmp( self, outdir, strbl, Nprocess ) : 
		freq = self.antarray.Ant.freq
		func, amp, color = ['FitBeam', 'FitVis'], [], ['b','r']
		try    : amp.append(self.Ampb.T)  # (Nv,Nf)
		except : amp.append(None)
		try    : amp.append(self.Ampv.T)
		except : amp.append(None)
		if (amp[0] is None and amp[1] is None) : 
			if (self.verbose) : print '    NOT exist  self.Ampb, self.Ampv'
			return
		if (self.verbose) : print '    Plotting    Amp     ......'
		send = tuple(amp) + (strbl,)
		bcast = (freq, func, color, outdir, 'amp')
		pool = jp.PoolFor(0, len(strbl), Nprocess)
		pool.map_async(_DoMultiprocess_Plot, send, bcast)



	def _PlotDeff( self, outdir, strbl, Nprocess, dyDeff ) : 
		freq = self.antarray.Ant.freq
		func, Deff, color = ['FitBeam', 'FitVis'], [], ['b','r']
		try    : Deff.append(300/(self.Sigmab.T*self.fwhm2sigmafactor*freq[None,:]))
		except : Deff.append(None)
		try    : Deff.append(300/(self.Sigmav.T*self.fwhm2sigmafactor*freq[None,:]))
		except : Deff.append(None)
		if (Deff[0] is None and Deff[1] is None and self.verbose) : 
			if (self.verbose) : print '    NOT exist  self.Sigmab, self.Sigmav'
			return
		if (self.verbose) : print '    Plotting    Deff    ......'
		send = tuple(Deff) + (strbl,)
		bcast = (freq, func, color, outdir, 'deff', dyDeff, self.fwhm2sigmafactor, self.antarray.Ant.dishdiam)
		pool = jp.PoolFor(0, len(strbl), Nprocess)
		pool.map_async(_DoMultiprocess_Plot, send, bcast)



	def _PlotLew( self, outdir, strbl, Nprocess, dyLew ) : 
		freq = self.antarray.Ant.freq
		func, Lew, color = ['FitPhase', 'FitVis'], [], ['b','r']
		try    : Lew.append(self.Lewp.T)  # (Nv,Nf)
		except : Lew.append(None)
		try    : Lew.append(self.Lewv.T)
		except : Lew.append(None)
		if (Lew[0] is None and Lew[1] is None) : 
			if (self.verbose) : print '    NOT exist  self.Lewp, self.Lewv'
			return
		if (self.verbose) : print '    Plotting    Lew     ......'
		bl =self.antarray.Blorder.baseline[self.antarray.visorder]
		send = tuple(Lew) + (strbl,)
		bcast = (freq, func, color, outdir, 'lew', bl, dyLew)
		pool = jp.PoolFor(0, len(strbl), Nprocess)
		pool.map_async(_DoMultiprocess_Plot, send, bcast)



	def _PlotPhaseadd( self, outdir, strbl, Nprocess ) : 
		freq = self.antarray.Ant.freq
	#	func,phaseadd,color = ['FitPhase','FitExp'], [], ['b','r']
		func, phaseadd, color = ['FitPhase','FitVis','FitExp'], [], ['b','r','c']
		try    : phaseadd.append(self.Phaseaddp.T*180/np.pi) # deg
		except : phaseadd.append(None)
		try    : phaseadd.append(self.Phaseaddv.T*180/np.pi)
		except : phaseadd.append(None)
		try    : phaseadd.append(self.Phaseadde.T*180/np.pi)
		except : phaseadd.append(None)
		if (phaseadd[0] is None and phaseadd[1] is None and phaseadd[2] is None) : 
			if (self.verbose) : print '    NOT exist  self.Phaseaddp, self.Phaseaddv, self.Phaseadde'
			return
		if (self.verbose) : print '    Plotting  Phaseadd  ......'
		bl =self.antarray.Blorder.baseline[self.antarray.visorder]
		send = tuple(phaseadd) + (strbl,)
		bcast = (freq, func, color, outdir, 'phaseadd', bl)
		pool = jp.PoolFor(0, len(strbl), Nprocess)
		pool.map_async(_DoMultiprocess_Plot, send, bcast)



	def _PlotFringe( self, array, label, color, lw=None, plotnv=None, plotfreq=None ) : 
		'''
		array: 
			vis to be plotted
			(1) ndarray
			(2) [ndarray, ndarray, ...], in this case, label, color, lw should has the same number

		label: 
			plt.plot(x, y, label=label)
			label must be str pair [label1, label2], 1=>real, 2=>imag

		color:
			color must be str pair [color1, color2], 1=>real, 2=>imag

		lw:
			linewidth, must be int, one int for real+imag
			can be lw=None
		'''
		# Plot, longest East-West baseline, center frequency
		if (self.timem[-1]-self.timem[0] <= 60) : xmajor = 5
		else : xmajor = 10
		#--------------------------------------------------
		istype = jp.IsType()
		if (not (istype.islist(array) or istype.istuple(array))) : 
			if (lw is None) : lw = 1 
			array, label, color, lw = [array], [label], [color], [lw]
		else : 
			if (lw is None) : lw = [1,1] 
		#--------------------------------------------------
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
		strbll = '=(%.3f, %.3f)' % tuple(self.antarray.Blorder.baseline[nvbl][:2])
		strbl = str(bl[0])+'-'+str(bl[1])
		nf = abs(freq-plotfreq)
		nf = np.where(nf==nf.min())[0][0]
		#--------------------------------------------------
		# Plot fringe
		vmax = []
		plt.figure(figsize=(17,6))
		for i in xrange(len(array)) : 
			vmax.append( abs(array[i][:,nf,nv]).max()*1.05 )
			scale = jp.SciNot(vmax[-1])[1]-2
			vmax[-1] /= 10.**scale
			plt.subplot(1,2,1)
			plt.plot(self.timem, array[i].real[:,nf,nv]/10.**scale, color=color[i][0], lw=lw[i], label=label[i][0])
			if (i == len(array)-1) : 
				plt.legend()
				plt.xlabel('time [min]', size=16)
				plt.xlim(self.timem.min(), self.timem.max())
				plt_axes('x', 'both', [xmajor,1])
				plt.ylabel('[A.U.]', size=16)
				vmax = np.array(vmax).max()
				plt.ylim(-vmax, vmax)
			plt.subplot(1,2,2)
			plt.plot(self.timem, array[i].imag[:,nf,nv]/10.**scale, color=color[i][1], lw=lw[i], label=label[i][1])
			if (i == len(array)-1) : 
				plt.legend()
				plt.xlabel('time [min]', size=16)
				plt.xlim(self.timem.min(), self.timem.max())
				plt_axes('x', 'both', [xmajor,1])
				plt.ylabel('[A.U.]', size=16)
				plt.ylim(-vmax, vmax)
		plt.suptitle('Fringe of baseline='+strbl+strbll+' @ '+freqstr+'MHz', fontsize=20)
		plt.savefig(self.outdir+'vis_'+strbl+'_'+freqstr+'MHz.png')
		plt.close()
		#--------------------------------------------------
		# Plot phase
		phase = np.angle(array[-1][:,nf,nv])
		phase = (phase*180/np.pi) %360
		plt.plot(self.timem, phase, 'bo', markersize=4)
		plt.xlabel('time [min]', size=16)
		plt.xlim(self.timem.min(), self.timem.max())
		plt_axes('x', 'both', [xmajor,1])
		plt.ylim(0, 360)
		plt.ylabel('[degree]', size=16)
		plt.title('Total phase of baseline='+strbl+strbll+' @ '+freqstr+'MHz', size=13)
		plt.savefig(self.outdir+'phase_'+strbl+'_'+freqstr+'MHz.png')
		plt.close()
