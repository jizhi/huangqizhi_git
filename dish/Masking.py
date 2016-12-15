import numpy as np
import multiprocessing
import jizhipy as jp
from jizhipy.Plot import *
import scipy.signal as spsn

from AntArray import *
##################################################




class NoiseSource( object ) : 


	def __init__( self, pixstart=None, pixlength=None, pixperiod=None ) : 
		''' NOTE THAT here set the pixel/index, not second ! 
		Generately, we set pixstart of the first file(filelist[0]). However, because inttime is not an integer, after long time, we may stager on pixel, in this case, we set pixstart of current (or one before current) Antarray (see self.Mask() below), then we may have a better result

		0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20
				. . . .                      .  .  .  .  
		. . .         . . . .  .  .  .  .              .  .

		pixstart=3,  pixlength=4 (3~6),  pixperiod=12 (15-3=12)
		'''
		self.pixstart, self.pixlength, self.pixperiod = int(round(pixstart)), int(round(pixlength)), int(round(pixperiod))


	def Mask( self, antarray ) : 
		'''
		(1) class:AntArray
		(2) np.ndarray with 3D (time, freq, baseline)
		'''
		N0 = antarray.Ant.N0
		nfile = antarray.File.nfile[0]
		Ntot = N0*nfile + antarray.vis.shape[0]
		#--------------------------------------------------
		n1 = np.arange(self.pixstart, Ntot, self.pixperiod)
		n1 = n1[n1>=N0*nfile]
		n2 = n1 + self.pixlength
		n2 = n2[n2<=Ntot]
		if (n2.size < n1.size) : n2 = np.append(n2, [Ntot])
		self.nmask = np.concatenate([n1[:,None], n2[:,None]], 1)-N0*nfile
		self.mask = np.zeros(antarray.vis.shape[0], bool)
		for i in xrange(len(self.nmask)) : 
			self.mask[self.nmask[i,0]:self.nmask[i,1]] = True
		return self.nmask, self.mask





##################################################
##################################################
##################################################





class Masking( object ) : 



	def __init__( self, antarray=None, verbose=True, outdir='' ) :
		'''
		antarray: 
			Must be instance of class:AntArray

		self.maskvalue:
			Resetd/good values
			Usage: antarray.vis[self.mask] = self.maskvalue to reset the bad values with good values
		'''
		try : 
			self.verbose = antarray.verbose
			self.antarray = antarray
		except : 
			self.verbose = bool(verbose)
			self.antarray = AntArray(verbose=False)
			self.antarray.verbose = self.verbose
		if (self.verbose) : print '-------------------- Masking --------------------\n'
		self.outdir = jp.Outdir((None,'file'), (0,'file'))
		try : self.outdir = jp.Mkdir(self.outdir + outdir)
		except : pass
		#--------------------------------------------------
		# Always have self.mask.shape == vis.shape
		try : 
			self.Nprocess = self.antarray.Nprocess
			self.mask = np.zeros(self.antarray.vis.shape, bool)
			self.maskvalue = np.array([], self.antarray.vis.dtype)
			self.outname = self.outdir + 'Masking_'+self.antarray.File.filedirname+'.hdf5'
		except : pass
		#	if (self.verbose) : jp.Raise(Warning, 'Wrong  antarray  in Masking.__init__()')
		self.masknoisesource, self.pixstart, self.pixlength, self.pixperiod = None, None, None, None





	def MaskNoiseSource( self, pixstart=None, pixlength=None, pixperiod=None ) :
		'''
		0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20
				. . . .                      .  .  .  .  
		. . .         . . . .  .  .  .  .              .  .

		pixstart=3,  pixlength=4 (3~6),  pixperiod=12 (15-3=12)

		self..pixstart
		self..pixlength
		self..pixperiod
		self.noisesource.nmask
		self.noisesource.mask
		self.masknoisesource  # same shape as self.mask

		(1) If can do pixstart+pixlength+pixperiod, then use them.
		(2) Else, calculate them automatically. In this case, times=pixstart, you can set it from 0 to 1 to make the result more correct.
		'''
		if (self.verbose) : 
			print 'Masking.MaskNoisesource'
		#	print '    start @', jp.Time(1)
		try : pixstart+pixlength+pixperiod
		except : 
			if (pixstart is None) : times = 0.5
			else : times = pixstart
			try : vis = abs(self.antarray.vis[:,self.antarray.vis.shape[1]/2,0].copy())
			except : vis = abs(self.antarray[:,self.antarray.shape[1]/2,0].copy())
			vis =np.arange(vis.size)[vis>=vis.max()*times] #@ times
			ngroup, nrange = jp.ArrayGroup(vis)
			n = nrange[:,1] - nrange[:,0]
			count = []
			while (n.size > 0) : 
				count.append( [n[n==n[0]].size, n[0]] )
				n = n[n!=n[0]]
			count = jp.Sort(np.array(count), '[:,0]', True)[0,1]
			n = nrange[:,1] - nrange[:,0]
			nrange = nrange[n==count]
			pixlength = nrange[0,1] - nrange[0,0]  # sure
			n = nrange[1:,0] - nrange[:-1,0]
			pixperiod = []
			while (n.size > 0) : 
				pixperiod.append( [n[n==n[0]].size, n[0]] )
				n = n[n!=n[0]]
			pixperiod = jp.Sort(np.array(pixperiod), '[:,0]', True)[0,1]
			n = nrange[1:,0] - nrange[:-1,0]
			n = np.where(n==pixperiod)[0][0]
			pixstart = nrange[n,0] + self.antarray.File.nfile[0]*self.antarray.Ant.N0
			pixstart = pixstart - pixstart/pixperiod * pixperiod
		if (self.verbose) : print '    pixstart=%i, pixlength=%i, pixperiod=%i' % (pixstart, pixlength, pixperiod)
		#--------------------------------------------------

		noisesource = NoiseSource(pixstart, pixlength, pixperiod)
		noisesource.Mask(self.antarray)
		self.masknoisesource = noisesource.mask[:,None,None] + np.zeros(self.mask.shape, bool)  # 3D
		self.pixstart, self.pixlength, self.pixperiod = noisesource.pixstart, noisesource.pixlength, noisesource.pixperiod
#		self.mask += self.masknoisesource
		''' self.mask is old, ==> self.maskvalue
		    self.masknoisesource is new, reset this part '''
		#--------------------------------------------------
		if (self.verbose) : 
			print '    Masked elements: %i in %i' % (self.masknoisesource.sum(), self.mask.size)
			print '    Reset the masked values ......'
		self._ResetMasked(self.masknoisesource)
	#	if (self.verbose) : print '    end   @', jp.Time(1)+'\n'
		if (self.verbose) : print





	def _ResetMasked( self, newmask ) : 
		'''
		Rest self.antarray.vis
		Update self.maskvalue, self.mask
		'''
		# self.maskvalue is good values
		maskback = self.mask.copy()
		badvalue = self.antarray.vis[maskback].copy()
		self.antarray.vis[self.mask] = self.maskvalue
		self.antarray.vis = np.ma.MaskedArray(self.antarray.vis, newmask)
		goodvalue = jp.ResetMasked(self.antarray.vis, 0, self.Nprocess)
		self.antarray.vis = self.antarray.vis.data
		self.antarray.vis[newmask] = goodvalue
		self.mask += newmask  #@#@
		self.maskvalue = self.antarray.vis[self.mask].copy()
		self.antarray.vis[maskback] = badvalue
		# self.maskvalue is bad values
	#	shape, size = self.mask.shape, self.mask.size
	#	vis = np.ma.MaskedArray(self.antarray.vis, newmask)
	#	newvalue = jp.ResetMasked(vis, 0, self.Nprocess)
	#	vis = self.antarray.vis
	#	oldvalue = vis[self.mask].copy()
	#	vis[self.mask] = self.maskvalue  # reset to original
	#	self.maskvalue = vis[self.mask + newmask]  #@#@
	#	vis[self.mask] = oldvalue
	#	vis[newmask] = newvalue
	#	self.mask += newmask  #@#@





	def MaskLoop( self, axis, per, times, nsigma=5, nloop=None, threshold=None, xhour=False ) : 
		'''
		(1) select every freq bins, mask along time
		(2) select every time bins, mask along freq
		(1) and (2) may get the same/similar self.mask array
		So, use (2) axis=0 (time) is enough !

		per, times:
			Use for smooth()
			In order to set a good per and times, you can try by hand and judge by eyes (plot and look at around the fringe/source)

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
		'''
		if (self.verbose) : 
			print 'Masking.MaskLoop'
			print '    start @', jp.Time(1)
		axis, per, times = int(round(axis)), int(round(per)), int(round(times))
		shape = self.antarray.vis.shape
		if (axis < 0) : axis = len(shape) + axis
		if (axis > len(shape)-1) : jp.Raise(Exception, 'axis='+str(axis)+' out of self.antarray.shape='+str(shape))
		if (per<=1 or times<=0) : return
		self.per, self.times = per, times
		#--------------------------------------------------

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
		if (self.verbose) : print ('    per=%i, times=%i, nsigma=%.1f, nloop=%i, threshold=%.3f' % (per, times, nsigma, nloop, threshold))
		#--------------------------------------------------

		# Move to axis=0
		if (axis != 0) : 
			vis = self.antarray.vis.copy()
			vis[self.mask] = self.maskvalue
			self.antarray.vis = jp.ArrayAxis(self.antarray.vis, axis, 0, 'move')
			vis = jp.ArrayAxis(vis, axis, 0, 'move')
			self.mask = jp.ArrayAxis(self.mask, axis, 0, 'move')
			self.maskvalue = vis[self.mask] 
			del vis
		#--------------------------------------------------

		if (self.verbose) : print '    done     before      after       diff     time'
		done, mask, Nmask = 0, 0, 0
		while (done < nloop) : 

			# Reset good values
			badvalue = self.antarray.vis[self.mask].copy()
			maskback = self.mask.copy()
			self.antarray.vis[maskback] = self.maskvalue
			#--------------------------------------------------

			before = ('%11i') % self.mask.sum()
			vis = self.antarray.vis
			dvis = vis - jp.Smooth(vis, 0, per, times, Nprocess=self.Nprocess)
			vistype = 'cross' if(dvis.dtype.name[:7]=='complex')else 'auto'
			#--------------------------------------------------

			if (vistype == 'cross') : 
				dvis.real, dvis.imag = dvis.real**2, dvis.imag**2
				sigma2 = np.ma.MaskedArray(dvis,self.mask).mean(0)
				newmask = (dvis.real>nsigma**2*sigma2.real) + (dvis.imag>nsigma**2*sigma2.imag) # mask
			#--------------------------------------------------

			elif (vistype == 'auto') : 
				dvis = dvis**2
				sigma2 = np.ma.MaskedArray(dvis, self.mask).mean(0)
				newmask = (dvis > nsigma**2*sigma2) # mask
			done += 1
			#--------------------------------------------------

			# ResetMasked
			self._ResetMasked(newmask)
			self.antarray.vis[maskback] = badvalue
			#--------------------------------------------------

			strdone = '%8i'    % done
			after   = ('%11i') % self.mask.sum()
			diff    = ('%11i') % (int(after)-int(before))
			if (self.verbose) : print strdone + before + after + diff +'    ', jp.Time(1)[11:]
			Nmask += newmask.sum()
			if (int(diff) <= dvis.size*threshold) : break
		#--------------------------------------------------

		if (self.verbose) : print '    Masked elements: %i in %i' % (Nmask, self.mask.size)
		self.PlotMask(None, None, xhour)
		if (axis != 0) : 
			vis = self.antarray.vis.copy()
			vis[self.mask] = self.maskvalue
			self.antarray.vis = jp.ArrayAxis(self.antarray.vis, 0, axis, 'move')
			vis = jp.ArrayAxis(vis, 0, axis, 'move')
			self.mask = jp.ArrayAxis(self.mask, 0, axis, 'move')
			self.maskvalue = vis[self.mask]
			del vis
		if (self.verbose) : print '    end   @', jp.Time(1)+'\n'





	def MaskManual( self, maskmanual ) : 
		if (maskmanual.shape != self.mask.shape) : jp.Raise(Exception, 'maskmanual.shape != self.mask.shape')
		if (self.verbose) : 
			print 'Masking.MaskManual'
			print '    start @', jp.Time(1)
			print '    Masked elements: %i in %i' % (self.maskmanual.sum(), self.mask.size)
			print '    Reset the masked values ......'
		self._ResetMasked(maskmanual)
		if (self.verbose) : print '    end   @', jp.Time(1)+'\n'





	def See( self, per=None, times=None, masked=False, masktime=None, outdir=None, show=False ) : 
		'''
		in See(), masktime must be pair: masktime=(n1, n2)
		masktime: just pixels, independent of self.inttime

		masked:
			True | False
			==True: This is MaskedArray, will plot original data and masked data
			==False: This is normal np.ndarray
		'''
		if (outdir is None) : outdir = self.outdir
		else : outdir = jp.Mkdir(outdir)  #@#@
	#	self.outdir = jp.Outdir((None,'file'), (0,'file'))
		try : per, times = int(round(per)), int(round(times))
		except : per, times = 0, 0
		nf, strfreq, nv, strblo, strbli, blname = self.antarray.Plotnfnv(None, None, False)
		if (self.verbose) : 
			print 'XXX.See'
			print '    per=%i, times=%i' % (per, times)
			print '    Plotting XXX.See_'+strfreq+'_'+strblo+'.png\n'
		#--------------------------------------------------

		if (not masked) : vis = self.antarray.vis[:,nf,nv].real
		else : 
			viso = self.antarray.vis[:,nf,nv].real
			vis = np.ma.MaskedArray(self.antarray.vis.copy(), self.mask)
			vis.data[self.mask] = self.maskvalue
			vis = vis[:,nf,nv].real.copy()  # masked
		#--------------------------------------------------
		if (per !=0 and times !=0) : 
			viss = jp.Smooth(vis.data, 0, per, times, Nprocess=self.Nprocess)
		#--------------------------------------------------

		if (masktime is not None) : 
			try : 
				n1, n2 = masktime
				if (vis[n1:n2].size == 0) : masktime = None
				else : plt.figure(figsize=(17, 6))
			except : masktime = None
		if (masktime is None) : plt.figure(figsize=(8,6))
		#--------------------------------------------------

		x = np.arange(vis.size)
		color = ['b', 'r']
		if (masked and (per!=0 and times!=0)) : color += ['g']
		if (masktime is not None) : plt.subplot(1,2,1)

		if (masked) : 
			plt.plot(x, viso, color=color[0], label='data.real')
			plt.plot(x, vis,  color=color[1], label='masked')
		else : 
			plt.plot(x, vis,  color=color[0], label='data.real')
		if (per != 0 and times != 0) : 
			plt.plot(x, viss, color=color[-1], label='smoothed\nper=%i, times=%i' % (per, times), lw=2)

		plt.legend(fontsize=12)
		plt.xlim(x[0], x[-1])
		plt.xlabel('Pixel', size=16)
		if (masktime is None) : plt.title('XXX.See, '+strfreq+', '+blname+strblo+strbli, size=14)
		#--------------------------------------------------

		if (masktime is not None) : 
			plt.subplot(1,2,2)
			if (masked) : 
				plt.plot(x[n1:n2], viso[n1:n2], color=color[0], label='data.real')
				plt.plot(x[n1:n2], vis[n1:n2],  color=color[1], label='masked')
			else : 
				plt.plot(x[n1:n2], vis[n1:n2],  color=color[0], label='data.real')
			if (per != 0 and times != 0) : 
				plt.plot(x[n1:n2], viss[n1:n2], color=color[-1], label='smoothed\nper=%i, times=%i' % (per, times), lw=2)

			plt.legend(fontsize=12)
			plt.xlim(x[n1], x[n2])
			plt.xlabel('Pixel', size=16)
			plt.title('masktime=(%i, %i)' % (n1, n2), size=16)
			plt.suptitle('XXX.See, '+strfreq+', '+blname+strblo+strbli, size=18)
		plt.savefig(outdir+'XXX.See_'+strfreq+'_'+strblo+'.png')
		if (show) : plt.show()
		else : plt.close()





	def PlotMask( self, nf, nv, xhour=False ) : 
		if (self.verbose and jp.SysFrame(0,2)[-3][-2]=='') : 
			print 'Masking.PlotMask'
		nf, strfreq, nv, strblo, strbli, blname = self.antarray.Plotnfnv(nf, nv)
		try : 
			if (xhour) : x = self.antarray.timem /60.
			else : x = self.antarray.timem
			xlabel = 'time [min]'
		except : 
			x = np.arange(self.antarray.vis.shape[0])
			xhour = False
			xlabel = 'Pixel'
		#--------------------------------------------------
		vistype = 'cross' if(self.antarray.vis.dtype.name[:7]=='complex')else 'auto'
		self.antarray.vis = np.ma.MaskedArray(self.antarray.vis, self.mask)

		for i in xrange(nf.size) : 
			for j in xrange(nv.size) : 
				plt.figure(figsize=(16, 12))
				vis = self.antarray.vis[:,nf[i],nv[j]]
				plt.subplot(2,1,1)
				plt.plot(x, vis.data, 'b-', label='data.real')
				plt.plot(x, vis, 'r-', label='RFI masked')
				plt.legend()
				plt.xlim(x.min(), x.max())
				if (xhour) : 
					plt_period(x, 24, 1)
					plt.xlabel("o'clock", size=16)
				else : plt.xlabel(xlabel, size=16)
				plt.title('Masking, '+strfreq[i]+', '+blname+strblo[j]+strbli[j], size=22)

				plt.subplot(2,1,2)  # zoom in
				plt.plot(x, vis.data, 'b-', label='data.real')
				plt.plot(x, vis, 'r-', label='RFI masked')
				plt.legend()
				plt.xlim(x.min(), x.max())
				if (xhour) : 
					plt_period(x, 24, 1)
					plt.xlabel("o'clock", size=16)
				else : plt.xlabel(xlabel, size=16)
			#	vmin1, vmax1 = vis.data.min(), vis.data.max()
				vmin2, vmax2 = vis.min(), vis.max()
				d = (vmax2 - vmin2)/4.
				vmin2, vmax2 = vmin2-d, vmax2+1*d
			#	vmin2, vmax2 = max(vmin2, vmin1), min(vmax2, vmax1)
				plt.ylim(vmin2, vmax2)
				plt.title('Zoom in', size=16)
				plt.savefig(self.outdir+'Masking_'+strfreq[i]+'_'+strblo[j]+'.png')
				plt.close()
				if (self.verbose) : print '    Plotting Masking_'+strfreq[i]+'_'+strblo[j]+'.png'

		self.antarray.vis = self.antarray.vis.data
		if (self.verbose and jp.SysFrame(0,2)[-3][-2]=='') : print





	def Save( self, outname=None ) : 
		''' outname: Absolute path of output .hdf5 '''
		if (outname is not None) : self.outname = str(outname)
		if (self.outname[-5:] !='.hdf5') : self.outname += '.hdf5'
		if (self.verbose) : 
			print 'Masking.Save'
			print '    Saving to  "'+jp.AbsPath(self.outname)+'"\n'
		masking = Masking(verbose=False)
		masking.__dict__ = self.__dict__
		masking.__dict__.pop('antarray')
		classhdf5 = jp.ClassHdf5(masking, self.outname, verbose=False)
		classhdf5.Save()





	def Read( self, outname=None ) : 
		''' outname: Absolute path of output .hdf5 '''
		if (outname is not None) : self.outname = str(outname)
		if (self.outname[-5:] !='.hdf5') : self.outname += '.hdf5'
		if (self.verbose) : 
			print 'Masking.Read'
			print '    Reading from  "'+jp.AbsPath(self.outname)+'"\n'
		masking = Masking(verbose=False)
		classhdf5 = jp.ClassHdf5(masking, self.outname, verbose=False)
		classhdf5.Read()
		masking.__dict__.pop('antarray') 
		self.__dict__.update(masking.__dict__)


