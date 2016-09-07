import numpy as np
import multiprocessing
import jizhipy as jp
from jizhipy.Plot import *
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


def _DoMultiprocess_vis( iterable ) : 
	vis, axis, per, times, nsigma = iterable
	print vis.shape, jp.Time(1)
	exit()
	if (per<=1 or times<=0) : return False
	dvis =vis.data -jp.Smooth(vis.data, axis, per, times)
	sigma = np.ma.MaskedArray(dvis**2, vis.mask).mean(axis)**0.5
	if (axis == 1) : mask = (abs(dvis) > nsigma*sigma[:,None,:])
	else : mask = (abs(dvis) > nsigma*sigma)
	return mask


class Masking( object ) : 

	def __init__( self, antarray=None ) : 
		if (antarray is None) : return
		self.antarray = antarray  #@ will it increast the memory?

		# Always have self.mask.shape == vis.shape
		shape = list(self.antarray.vis.shape)
		shape[2] = len(self.antarray.visorder)
		self.mask = np.zeros(shape, bool)
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
		self.noisesource = NoiseSource(pixstart, pixlength, pixperiod)
		self.noisesource.Mask(self.antarray)
		self.masknoisesource = self.noisesource.mask[:,None,None]
		self.mask += self.masknoisesource



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



	def MaskLoop( self, timeper=60, timetimes=1, freqper=3, freqtimes=1, nsigma=5, nloop=None, threshold=None, multipool=True ) : 
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

		if nloop!=None, use nloop and don't use threshold
		elif nloop==None but threshold!=None, use threshold
		elif nloop==None and threshold==None, stop with threshold=0.001

		multipool:
			True/False

		return self.maskloop
		Total mask: self.maskloop + self.masknoisesource
		'''
		# Read the whole real data as MaskedArray
		print 'antarray.vis[:,:,visorder]', jp.Time(1)
		vis = self.antarray.vis[:,:,self.antarray.visorder]
		if (self.antarray.visorder.size==1) : vis = vis[:,:,None]
		vis = np.ma.MaskedArray(vis, self.mask)
		print 'antarray.vis[:,:,visorder]  -->  END', jp.Time(1)
		# Convenient to Smooth(), using non-masked array is much faster than masked array
	#	vis = vis[:16600, 307-20:307+20]  # PAON4 test
	#	vis = vis[:,:3]  # PAON4 test
		vis.data[vis.mask] = vis.mean() # data

		x = np.arange(vis.shape[0])
		plt.plot(x, vis.real[:,1,0], 'b-')

		if (nloop) : 
			stop, threshold = 1, 0
			if (nloop <= 0) : nloop = 1
		else : 
			stop, nloop = 2, 100
			if (not threshold) : threshold = 0.001

		multipool = False
		done, mask = 0, 0
		while (done < nloop) : 
			if (done == 0) : 
				mask = np.concatenate([vis.mask, vis.mask], 2)
				vis=np.concatenate([vis.data.real, vis.data.imag],2)
				vis = np.ma.MaskedArray(vis, mask)
				Nv = vis.shape[-1]
			else : vis.mask = np.concatenate([mask, mask], 2)
			masksum, mask = vis.mask.sum()/2, 0

			if (not multipool) : 
				#---------- Single thread START ----------
				# time
				if (timeper>=2 and timetimes>=1) : 
					print '---1'
					dvis = vis.data-jp.Smooth(vis.data, 0, timeper, timetimes, Nprocess=None)
					print '---2'
					sigma = np.ma.MaskedArray(dvis**2, vis.mask).mean(0)**0.5
					print '---3'
					vis.mask += (abs(dvis) > nsigma*sigma) # mask
					print '---4'
					print jp.Time(1)
				# freq
				if (freqper>=2 and freqtimes>=1) : 
					dvis = vis.data - jp.Smooth(vis.data, 1, freqper, freqtimes, Nprocess=None)
					print '---5'
					sigma = np.ma.MaskedArray(dvis**2, vis.mask).mean(1)**0.5
					print '---6'
					vis.mask += (abs(dvis) > nsigma*sigma[:,None,:])
					print '---7'
					print jp.Time(1)
				#---------- Single thread  END  ----------

			else : 
				#---------- Pool START ----------
				# Can't pass very large array (vis) to other progress
				# Pool method is not good!
				pool = multiprocessing.Pool(2)
				mask = pool.map_async(_DoMultiprocess_vis, ((vis, 0, timeper, timetimes, nsigma), (vis, 1, freqper, freqtimes, nsigma))).get(10**10)
				vis.mask += mask[0] + mask[1]
				#---------- Pool  END  ----------


		#	plt.plot(x, vis.real[:,1,0], 'r-')
		#	plt.show()
		#	jp.Raise()

			mask = vis.mask[:,:,:Nv/2] + vis.mask[:,:,Nv/2:]
			done += 1
			diff = mask.sum()-masksum

			print done, '  ', masksum, '  ', mask.sum(), '  ', diff, '  ', jp.Time(1)
			jp.Raise()
			if (stop==2 and diff<=vis.size/2*threshold) : break

		try : self.maskloop = mask - self.masknoisesource
		except : self.maskloop = mask
		self.mask += mask



	def MaskManual( self, maskmanual ) : 
		if (maskmanual.shape != self.mask.shape) : jp.Raise(Exception, 'maskmanual.shape != self.mask.shape')
		self.maskmanual = maskmanual
		self.mask += maskmanual


##################################################
##################################################
##################################################

