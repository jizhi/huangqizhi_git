import numpy as np
import multiprocessing
import jizhipy as jp
from jizhipy.Plot import *
import scipy.signal as spsn
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


	def __init__( self, antarray=None, Nprocess=None, verbose=True ) : 
		self.Nprocess = Nprocess
		self.verbose = verbose
		if (antarray is None) : return
		self.antarray = antarray  #@ will it increast the memory?
		# Always have self.mask.shape == vis.shape
		shape = list(self.antarray.vis.shape)
		shape[2] = len(self.antarray.visorder)
		self.mask = np.zeros(shape, bool)
		#----------
		vistype = self.antarray.vistype[:-1]
		value = np.array(1, self.antarray.vis.dtype)
		if (vistype == 'auto') : dtype = value.real.dtype
		else : dtype = value.dtype
		self.maskvalue = np.array([], dtype)
		self.outdir = jp.Outdir((None,'file'), (0,'file'))



	def MaskNoiseSource( self, pixstart, pixlength, pixperiod ) :
		'''
		self.noisesource.pixstart
		self.noisesource.pixlength
		self.noisesource.pixperiod
		self.noisesource.nmask
		self.noisesource.mask
		self.masknoisesource
		'''
		if (self.verbose) : print 'Masking.MaskNoisesource: start @', jp.Time(1)
		self.noisesource = NoiseSource(pixstart, pixlength, pixperiod)
		self.noisesource.Mask(self.antarray)
		self.masknoisesource = self.noisesource.mask[:,None,None] + np.zeros(self.mask.shape, bool)  # 3D
		#--------------------------------------------------
		# Read and filter vis
		vistype = self.antarray.vistype[:-1]
		vis = self.antarray.vis[:,:,self.antarray.visorder]
		if (vistype == 'auto') : vis = vis.real
		if (len(vis.shape) == 2) : vis = vis[:,:,None]
		shape = vis.shape
		#--------------------------------------------------
		maskidx = np.arange(vis.size).reshape(vis.shape)[self.mask]  # ==> self.maskvalue
		self.mask += self.masknoisesource
		maskother = self.mask - self.masknoisesource
		#--------------------------------------------------
		# other
		maskotheridx = np.arange(vis.size).reshape(vis.shape)[maskother]
		maskothervalue = vis.flatten()*0
		maskothervalue[maskidx] = self.maskvalue
		maskothervalue = maskothervalue[maskotheridx]
		# noisesource
		masknoisesourceidx = np.arange(vis.size).reshape(vis.shape)[self.masknoisesource]
		# Reset
		vis = vis.flatten()
		vis[maskotheridx] = maskothervalue
		vis = np.ma.MaskedArray(vis.reshape(shape), self.masknoisesource)
		vis = jp.ResetMasked(vis, 0, self.Nprocess)
		masknoisesourcevalue = vis[self.masknoisesource]
		# Total
		maskidx = np.concatenate([maskotheridx, masknoisesourceidx])
		maskvalue = np.concatenate([maskothervalue, masknoisesourcevalue])
		maskvalue = np.array([maskidx, maskvalue], self.maskvalue.dtype)
		self.maskvalue = jp.Sort(maskvalue, '[0,:]')[1]
		if (self.verbose) : print 'Masking.MaskNoisesource:  end  @', jp.Time(1)+'\n'



	def See( self, freqpix, timelim, timeper, timetimes, timepix, freqlim, freqper, freqtimes, show=False ) : 
		'''
		Select the longest East-West baseline
		(1) See the fringe
		(2) See the frequency response

		freqpix, timepix:
			Can be int: use this one
			Or list/tuple/ndarray of int: average over these

		timelim, freqlim:
			int pair: (int,int) or [int,int] or np.array([int,int])
			Range of the figure
			==None means full range
		'''
		# freqpix, timepix
		freqpix, timepix = np.array(freqpix).round().astype(int), np.array(timepix).round().astype(int)
		strfreqpix = str(freqpix.min())
		if (freqpix.size!=1) : strfreqpix += '-'+str(freqpix.max())
		strtimepix = str(timepix.min())
		if (timepix.size!=1) : strtimepix += '-'+str(timepix.max())
		# timelim, freqlim
		try : 
			timelim1, timelim2 = np.array(timelim)[:2].round().astype(int)
			strtimelim = str(timelim1)+'-'+str(timelim2)
		except : timelim, strtimelim = None, 'None'
		try : 
			freqlim1, freqlim2 = np.array(freqlim)[:2].round().astype(int)
			strfreqlim = str(freqlim1)+'-'+str(freqlim2)
		except : freqlim, strfreqlim = None, 'None'
		# per, times
		timeper = np.array(timeper).round().astype(int).min()
		timetimes = np.array(timetimes).round().astype(int).min()
		freqper = np.array(freqper).round().astype(int).min()
		freqtimes = np.array(freqtimes).round().astype(int).min()
		# t-Amp
		bl = abs(self.antarray.Blorder.baseline[self.antarray.visorder][:,0])
		bl = np.where(bl==bl.max())[0][0]
		nbl = self.antarray.visorder[bl]
		strbl = '(%.3f, %.3f, %.3f)' % tuple(self.antarray.Blorder.baseline[nbl])
		vist = self.antarray.vis[:,freqpix,nbl].real
		if (len(vist.shape) != 1) : vist = vist.mean(1)
		# f-Amp
		na = self.antarray.Blorder.blorder[nbl][0]
		na = self.antarray.Blorder.Bl2Order(na, na)
		visf = self.antarray.vis[timepix,:,na].real
		if (len(visf.shape) != 1) : visf = visf.mean(0)
		t = np.arange(vist.size)
		f = np.arange(visf.size)
		plt.figure(figsize=(17,6))
		plt.subplot(1,2,1)
		plt.plot(t, vist, 'b-', label='t-Amp')
		plt.plot(t, jp.Smooth(vist, 0, timeper, timetimes), 'r-', label='timeper='+str(timeper)+', timetimes='+str(timetimes))
		plt.legend()
		if (timelim) : plt.xlim(timelim[0], timelim[1])
		else : plt.xlim(t.min(), t.max())
		plt.xlabel('t-points', size=16)
		plt.title('t-Amp', size=16)
		plt.subplot(1,2,2)
		plt.plot(f, visf, 'b-', label='f-Amp')
		plt.plot(f, jp.Smooth(visf, 0, freqper, freqtimes), 'r-', label='freqper='+str(freqper)+', freqtimes='+str(freqtimes))
		plt.legend()
		if (freqlim) : plt.xlim(freqlim[0], freqlim[1])
		else : plt.xlim(f.min(), f.max())
		plt.xlabel('f-points', size=16)
		plt.title('f-Amp', size=16)
		plt.suptitle('nbl='+str(nbl)+', bl='+strbl, size=16)
		figname = 'Masking.See_'+strfreqpix+'.'+strtimelim+'.'+str(timeper)+'.'+str(timetimes)+'_'+strtimepix+'.'+strfreqlim+'.'+str(freqper)+'.'+str(freqtimes)+'.png'
		plt.savefig(self.outdir+figname)
		if (show) : plt.show()
		plt.close()



	def MaskLoop( self, timeper=60, timetimes=1, nsigma=5, nloop=None, threshold=None, filtersize=None ) : 
		'''
		timeper, timetimes, freqper, freqtimes:
			Use for smooth()
			In order to set a good per+times, you can try by hand and judge by eyes (plot and look at around the fringe/source)

		nsigma:
			Calculate sigma of each freq and vis along time (sigma.shape=(Nf, Nv)), value > nsigma*sigma will be considered as RFI

		nloop:
			Each time we will mask the value>nsigma*sigma, and then recalculate a new sigma, and do again. nloop set how many times we will do.
			If set nloop=None, stop with threshold

		threshold:
			masknew.sum()-maskold.sum() < vis.size*threshold, stop

		If nloop=None, use threshold
		If threshold=None, use nloop
		If nloop!=None and threshold!=None, use one that satisfies first
		If nloop==None and threshold==None, set nloop=10, threshold=0.001

		filtersize:
			A scalar: filtersize=5 
			List of int for time-axis and freq-axis: [61, 5]
			If ==None: filtersize=[timeper, freqper]

		return self.maskloop
		Total mask: self.maskloop + self.masknoisesource
		'''
		if (self.verbose) : print 'Masking.MaskLoop: start @', jp.Time(1)
		timeper, timetimes = np.array([timeper, timetimes]).round().astype(int)
		if (timeper<=1 or timetimes<=0) : return
		if (nloop) : 
			try : nloop = int(round(nloop))
			except : nloop = 1
			if (nloop <= 0) : nloop = 1
		if (threshold) : 
			try : threshold+0
			except : threshold = 0.001
			if (threshold < 0) : threshold = 0
		if (nloop and not threshold) : threshold, strthreshold = 0, 'None'
		elif (not nloop and threshold) : nloop, strnloop = 100, 'None'
		elif (not nloop and not threshold) : nloop, threshold, strnloop, strthreshold = 10, 0.001, 'None', 'None'
		else : nloop, threshold
		#--------------------------------------------------
		# Read data
		vistype = self.antarray.vistype[:-1]
		vis = self.antarray.vis[:,:,self.antarray.visorder]
		# 3D, auto-real, cross-complex
		if (vistype == 'auto') : vis = vis.real
		if (len(vis.shape) == 2) : vis = vis[:,:,None]
		shape = vis.shape
		# Convenient to Smooth(), using non-masked array is much faster than masked array
		# ntimefringe=0-16600, nfreq1400=307 for PAON4
		vis[self.mask] = self.maskvalue  #@ Reset
		vis = np.ma.MaskedArray(vis, self.mask)
		#--------------------------------------------------
		maskidx=[np.arange(vis.size).reshape(vis.shape)[self.mask]]
		maskvalue = [self.maskvalue]
		if (self.verbose) : print '    done     before      after       diff     time'
		done, mask = 0, 0
		while (done < nloop) : 
			maskback = vis.mask.copy()
			#--------------------
			dvis = vis.data - jp.Smooth(vis.data, 0, timeper, timetimes, Nprocess=self.Nprocess)
			if (vistype == 'cross') : 
				dvis.real, dvis.imag=dvis.real**2, dvis.imag**2
				sigma2=np.ma.MaskedArray(dvis,vis.mask).mean(0)
				vis.mask += (dvis.real>nsigma**2*sigma2.real) + (dvis.imag>nsigma**2*sigma2.imag) # mask
			elif (vistype == 'auto') : 
				dvis = dvis**2
				sigma2=np.ma.MaskedArray(dvis,vis.mask).mean(0)
				vis.mask += (dvis > nsigma**2*sigma2) # mask
			done += 1
			# self.maskvalue
			masknew = vis.mask - maskback
			visreset = np.ma.MaskedArray(vis.data, masknew)
			visreset = jp.ResetMasked(visreset, 0, self.Nprocess)
			maskidx.append(np.arange(vis.size).reshape(vis.shape)[masknew])
			maskvalue.append(visreset[masknew].copy())
			visreset = 0 #@
			# Reset masked
			vis = vis.flatten()
			vis.data[maskidx[-1]] = maskvalue[-1]
			vis = vis.reshape(shape)
			if (self.verbose) : 
				strdone = '%8i'    % done
				before  = ('%11i') % maskback.sum()
				after   = ('%11i') % vis.mask.sum()
				diff    = ('%11i') % (vis.mask.sum()-maskback.sum())
				print strdone + before + after + diff +'    ', jp.Time(1)[11:]
			if (vis.mask.sum()-maskback.sum() <= vis.size*threshold) : break
		#--------------------------------------------------
		maskidx = np.concatenate(maskidx)
		maskvalue = np.concatenate(maskvalue)
		self.maskvalue = jp.Sort(np.array([maskidx, maskvalue]), '[0,:]')[1]
		self.mask = vis.mask
		try    : self.maskloop = self.mask - self.masknoisesource
		except : self.maskloop = self.mask
		if (self.verbose) : print 'Masking.MaskLoop:  end  @', jp.Time(1)+'\n'



	def MaskManual( self, maskmanual ) : 
		if (maskmanual.shape != self.mask.shape) : jp.Raise(Exception, 'maskmanual.shape != self.mask.shape')
		if (self.verbose) : print 'Masking.MaskManual: start @', jp.Time(1)
		maskback = self.mask.copy()
		self.mask += maskmanual
		self.maskmanual = self.mask - maskback
		#--------------------------------------------------
		maskidx = [np.arange(vis.size).reshape(vis.shape)[maskback]]
		maskvalue = [self.maskvalue]
		maskidx.append(np.arange(self.mask.size).reshape(self.mask.shape)[self.maskmanual])
		#--------------------------------------------------
		vistype = self.antarray.vistype[:-1]
		vis = self.antarray.vis[:,:,self.antarray.visorder]
		if (vistype == 'auto') : vis = vis.real
		if (len(vis.shape) == 2) : vis = vis[:,:,None]
		vis[maskback] = self.maskvalue
		vis = np.ma.MaskedArray(vis, self.maskmanual)
		vis = jp.ResetMasked(vis, 0, self.Nprocess)
		maskvalue.append(vis[self.maskmanual])
		#--------------------------------------------------
		maskvalue = np.array([np.concatenate(maskidx), np.concatenate(maskvalue)])
		self.maskvalue = jp.Sort(maskvalue, '[0,:]')[1]
		if (self.verbose) : print 'Masking.MaskManual:  end  @', jp.Time(1)+'\n'



##################################################
##################################################
##################################################


