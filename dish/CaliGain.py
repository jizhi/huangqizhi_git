import os
import numpy as np
import jizhipy as jp
from jizhipy.Plot import *
import scipy.signal as spsn
##################################################



class CaliGain( object ) : 


	def __init__(self, antarray=None, masking=None, verbose=True):
		self.verbose = verbose
		if (antarray is not None) : self.antarray = antarray
		if (masking  is not None) : self.masking  = masking
		self.outdir = jp.Outdir((None,'file'), (0,'file'))


	def Window( self, size ) : 
		''' How many time points to calculate one std? '''
		N = int(round(self.antarray.vis.shape[0] / size))
		self.nsplit = np.linspace(0, self.antarray.vis.shape[0], N+1).round().astype(int)



	def Gainnu( self ) : 
		''' Calculate g(nu) and plot '''
		outdir = self.outdir + 'Gainnu/'
		if (not os.path.exists(outdir)) : os.makedirs(outdir)
		if (self.verbose) : 
			print 'CaliGain.Gainnu: start @', jp.Time(1)
			progressbar = jp.ProgressBar('    completed:', len(self.nsplit)-1+self.antarray.visorder.size)
		# Read data
		vistype = self.antarray.vistype[:-1]
		if (vistype != 'auto') : 
			print '    Error: CaliGaint.Gainnu() should use "auto-correlation", but now AntArray.vistype="'+self.antarray.vistype+'"'
			if (self.verbose) : print 'CaliGain.Gainnu:  end  @', jp.Time(1)+'\n'
			exit()
		#--------------------------------------------------
		vis = self.antarray.vis[:,:,self.antarray.visorder].real
		if (len(vis.shape) == 2) : vis = vis[:,:,None]  # 3D
		# Replace masks
		vis[self.masking.mask] = self.masking.maskvalue
		# Calculate
		self.gainnu = np.zeros((self.nsplit.size-1,)+vis.shape[1:], vis.dtype)
		for i in xrange(len(self.nsplit)-1) : 
			if (self.verbose) : progressbar.Progress()
			self.gainnu[i] = vis[self.nsplit[i]:self.nsplit[i+1]].mean(0)
		#--------------------------------------------------
		# Plot
		gnusmooth = jp.Smooth(self.gainnu, 1, self.gainnu.shape[1]/10, 1)
		dgnu = (self.gainnu - gnusmooth).std(1)
		freq = self.antarray.Ant.freq
		color = plt_color(self.gainnu.shape[0])
		for i in xrange(self.gainnu.shape[2]) : 
			if (self.verbose) : progressbar.Progress()
			chanidx = self.antarray.Blorder.blorder[self.antarray.visorder][i][0]
			j = np.where(dgnu[:,i]==dgnu[:,i].min())[0][0]
			vmin, vmax = gnusmooth[j,:,i].min(), gnusmooth[j,:,i].max()
			vmin, vmax = vmin-0.05*(vmax-vmin), vmax+0.2*(vmax-vmin)
			if (vmin < 0) : vmin = 0
			for j in xrange(self.gainnu.shape[0]) : 
				plt.plot(freq, self.gainnu[j,:,i], color=color[j], label=str(j+1))
		#	plt.legend(fontsize=10)
			plt.xlabel(r'$\nu$ [MHz]', size=16)
			plt.xlim(int(round(freq.min())), int(round(freq.max())))
			plt_axes('x', 'both', [25,5])
			plt.ylabel('[A.U.]', size=16)
			plt.ylim(vmin, vmax)
			plt.title(r'$g(\nu)$ of channel '+str(chanidx)+', total '+str(len(color)+1)+' curves', size=16)
			plt.savefig('haha.png')
			plt.savefig(outdir+'gnu_'+str(chanidx)+'.png')
			plt.close()
		if (self.verbose) : print 'CaliGain.Gainnu:  end  @', jp.Time(1)+'\n'







	def _nfreq( self, nfreq=None ) : 
		'''
		nfreq:
			(1) int: freq bin
			(2) str with 'Hz': '1420.406MHz' or '1.420406GHz' or '1420406000Hz', etc
			(3) None: center bin

		return [nfreq, strfreq]
			nfreq: AntArray.Ant.freq[nfreq], AntArray.vis[:,nfreq]
		'''
		istype = jp.IsType()
		try : 
			nfreq = int(round(nfreq))
			strfreq = '%.1fMHz' % self.antarray.Ant.freq[nfreq]
		except : 
			if (istype.isstr(nfreq) and 'hz' in nfreq.lower()) : 
				strfreq = nfreq
				n = nfreq.lower().rfind('hz')
				if (nfreq[n-1].lower() == 'm') : nfreq = float(nfreq[:n-1])  # MHz
				elif (nfreq[n-1].lower() == 'g') : nfreq = float(nfreq[:n-1])*1000
				else : nfreq = float(nfreq[:n-1])*1e-6
				nfreq = abs(self.antarray.Ant.freq-nfreq)
				nfreq = np.where(nfreq==nfreq.min())[0][0]
			else : 
				nfreq = self.antarray.vis.shape[1] /2
				strfreq = '%.1fMHz' % self.antarray.Ant.freq[nfreq]
		return [nfreq, strfreq]


	def _nbl( self, nbl=None ) : 
		'''
		nbl:
			(1) int: order of bl
			(2) int pair (int1, int2): channel pair of baseline
			(3) None: longest East-West baseline
		return [nbl, strbl, norder]

			nbl: AntArray.Blorder.blorder[nbl], AntArray.vis[:,:,nbl]
			norder: AntArray.visorder[norder], AntArray.vis[:,:,AntArray.visorder[norder]]
		'''
		vistype = self.antarray.vistype[:-1]
		try : 
			nbl = int(round(nbl))
			if (nbl not in self.antarray.visorder) : raise
		except : 
			try : 
				nbl = self.antarray.Blorder.Bl2Order(nbl[:2])
				if (nbl not in self.antarray.visorder) : raise
			except : 
				bl = abs(self.antarray.Blorder.baseline[self.antarray.visorder][:,0])
				bl = np.where(bl==bl.max())[0][0]
				nbl = self.antarray.visorder[bl]
		norder = self.antarray.Blorder.blorder[nbl]
		bl     = self.antarray.Blorder.baseline[nbl]
		if (vistype == 'auto') : strbl = 'auto:%i' % norder[0]
		else : strbl = 'cross:%i-%i=(%.3f, %.3f, %.3f)' % (tuple(norder)+tuple(bl))
		norder = np.where(self.antarray.visorder==nbl)[0][0]
		return [nbl, strbl, norder]



	def See( self, timeper, timetimes, nfreq=None, nbl=None ) : 
		'''
		nfreq:
			(1) int: freq bin
			(2) str with 'Hz': '1420.406MHz' or '1.420406GHz' or '1420406000Hz', etc
			(3) None: center bin

		nbl:
			(1) int: order of bl
			(2) int pair (int1, int2): channel pair of baseline
			(3) None: longest East-West baseline
		'''
		if (self.verbose) : print 'CaliGain.See: start @', jp.Time(1)
		vistype = self.antarray.vistype[:-1]
		nbl, strbl, norder = self._nbl(nbl)
		nfreq, strfreq = self._nfreq(nfreq)
		# time
		timeper, timetimes = np.array([timeper, timetimes]).round().astype(int)
		# Read data
		vis = self.antarray.vis[:,:,self.antarray.visorder]
		if (len(vis.shape) == 2) : vis = vis[:,:,None]
		if (vistype == 'auto') : vis = vis.real
		# Replace masks
		vis[self.masking.mask] = self.masking.maskvalue
		# Select
		vis = vis[:,nfreq,norder]
		# Get the fluctuation
		viss = jp.Smooth(vis, 0, timeper, timetimes)
		# Plot
		x = np.arange(vis.size)
		if (vistype == 'auto') : 
			plt.figure(figsize=(8,8))
			plt.subplot(2,1,1)
			plt.plot(x, vis, 'b-', label='data')
			plt.plot(x, viss, 'r-', label='smoothed')
			plt.legend()
			plt.xlabel('time', size=16)
			plt.xlim(x.min(), x.max())
			plt.ylabel('[A.U.]', size=16)
			plt.subplot(2,1,2)
			plt.plot(x, vis-viss, 'g-', label='residual')
			plt.legend()
			plt.xlabel('time', size=16)
			plt.xlim(x.min(), x.max())
			plt.ylabel('[A.U.]', size=16)
		else : 
			plt.figure(figsize=(17,8))
			plt.subplot(2,2,1)
			plt.plot(x, vis.real, 'b-', label='data.real')
			plt.plot(x, viss.real, 'r-', label='smoothed.real')
			plt.legend()
			plt.xlabel('time', size=16)
			plt.xlim(x.min(), x.max())
			plt.ylabel('[A.U.]', size=16)
			plt.subplot(2,2,3)
			plt.plot(x, vis.real-viss.real, 'g-', label='residual.real')
			plt.legend()
			plt.xlabel('time', size=16)
			plt.xlim(x.min(), x.max())
			plt.ylabel('[A.U.]', size=16)
			plt.subplot(2,2,2)
			plt.plot(x, vis.imag, 'b-', label='data.imag')
			plt.plot(x, viss.imag, 'r-', label='smoothed.imag')
			plt.legend()
			plt.xlabel('time', size=16)
			plt.xlim(x.min(), x.max())
			plt.ylabel('[A.U.]', size=16)
			plt.subplot(2,2,4)
			plt.plot(x, vis.imag-viss.imag, 'g-', label='residual.imag')
			plt.legend()
			plt.xlabel('time', size=16)
			plt.xlim(x.min(), x.max())
			plt.ylabel('[A.U.]', size=16)
		plt.suptitle('Flucation, '+strfreq+', '+strbl, size=20)
		plt.savefig(self.outdir+'CaliGain.See_'+strfreq+'_'+strbl.split('=')[0]+'.png')
		plt.close()
		if (self.verbose): print 'CaliGain.See:  end  @', jp.Time(1)+'\n'



	def Gaint(self, timeper=3,timetimes=100, nfreq=None,nbl=None):
		'''
		Actually the return is:
			sigma_noise * gain(t)

		nfreq, nbl:
			Will plot a figure with nfreq, nbl

		nfreq:
			(1) int: freq bin
			(2) str with 'Hz': '1420.406MHz' or '1.420406GHz' or '1420406000Hz', etc
			(3) None: center bin

		nbl:
			(1) int: order of bl
			(2) int pair (int1, int2): channel pair of baseline
			(3) None: longest East-West baseline
		'''
		if (self.verbose) : 
			print 'CaliGain.Gaint: start @', jp.Time(1)
			progressbar = jp.ProgressBar('completed:', len(self.nsplit)-1)
		# Read data
		vistype = self.antarray.vistype[:-1]
		vis = self.antarray.vis[:,:,self.antarray.visorder]
		if (vistype == 'auto') : vis = vis.real
		if (len(vis.shape) == 2) : vis = vis[:,:,None]  # 3D
		# Replace masks
		vis[self.masking.mask] = self.masking.maskvalue
		# Get the fluctuation
		vis -= jp.Smooth(vis, 0, timeper, timetimes)
	#	vis = np.ma.MaskedArray(vis, self.masking.mask)
		# Calculate the std of noise
		self.gaint = np.zeros((len(self.nsplit)-1,)+vis.shape[1:], vis.dtype)
		for i in xrange(len(self.nsplit)-1) : 
			if (self.verbose) : progressbar.Progress()
			if (vistype == 'auto') : self.gaint[i] = vis[self.nsplit[i]:self.nsplit[i+1]].std(0)
			else : 
				self.gaint[i].real = vis.real[self.nsplit[i]:self.nsplit[i+1]].std(0)
				self.gaint[i].imag = vis.imag[self.nsplit[i]:self.nsplit[i+1]].std(0)
		#---------- Plot ----------
		nbl, strbl, norder = self._nbl(nbl)
		nfreq, strfreq = self._nfreq(nfreq)
		gaint = self.gaint[:,nfreq,norder]
		# ylim
		gaintm = spsn.medfilt(gaint.real, 7)
		ymax, ymin = gaintm.max()*1.1, gaintm.min()*0.9
		self._plot = (nbl,strbl,norder, nfreq,strfreq, ymax,ymin)
		# Plot
		x = np.arange(gaint.size)
		if (vistype == 'auto') : 
			plt.plot(x, gaint, 'bo', markersize=5, label='auto')
			plt.legend()
			plt.xlim(x.min(), x.max())
			plt.ylim(ymin, ymax)
			plt.xlabel('time', size=16)
			plt.ylabel('[A.U.]', size=16)
		else : 
			plt.plot(x, gaint.real, 'bo', markersize=5,label='real')
			plt.plot(x, gaint.imag, 'ro', markersize=5,label='imag')
			plt.legend()
			plt.xlim(x.min(), x.max())
			plt.ylim(ymin, ymax)
			plt.xlabel('time', size=16)
			plt.ylabel('[A.U.]', size=16)
		plt.title(r'$G(t) \cdot \sigma_n$, '+strfreq+', '+strbl, size=16)
		plt.savefig(self.outdir+'CaliGain.Gaint_'+strfreq+'_'+strbl.split('=')[0]+'.png')
		if (self.verbose) : print 'CaliGain.Gaint:  end  @', jp.Time(1)+'\n'



	def Smooth( self, filtersize, smoothtimes ) : 
		'''
		Smooth self.gaint

		filtersize:
			int: use this int for both time-axis and freq-axis
			int pair: (int1, int2), int1 for time-axis, int2 for freq-axis

		smoothtimes:
			int: for jp.Smooth() over time-axis
		'''
		if (self.verbose) : print 'CaliGain.Smooth: start @', jp.Time(1)
		# filtersize
		try : filtersize = filtersize[:2]
		except : filtersize = (5, 3)
		if (len(filtersize) == 1) : filtersize = (filtersize[0], filtersize[0])
		# smoothtimes
		try : smoothtimes = int(round(smoothtimes))
		except : smoothtimes = self.gaint.shape[0]/5
		# Do smoothing
		vistype = self.antarray.vistype[:-1]
		self.gaintsmooth = self.gaint*0
		for i in xrange(self.gaint.shape[2]) : 
			if (vistype == 'auto') : self.gaintsmooth[:,:,i] = spsn.medfilt2d(self.gaint[:,:,i], filtersize)
			else : 
				self.gaintsmooth[:,:,i].real = spsn.medfilt2d(self.gaint[:,:,i].real, filtersize)
				self.gaintsmooth[:,:,i].imag = spsn.medfilt2d(self.gaint[:,:,i].imag, filtersize)
		self.gaintsmooth = jp.Smooth(self.gaint, 0, 3, smoothtimes)
		#--------------------------------------------------
		#---------- Plot ----------
		if ('_plot' in self.__dict__.keys()) : 
			nbl,strbl,norder, nfreq,strfreq, ymax,ymin = self._plot
			gaint  = self.gaint[:,nfreq,norder]
			gaints = self.gaintsmooth[:,nfreq,norder]
			# Plot
			x = np.arange(gaint.size)
			if (vistype == 'auto') : 
				plt.plot(x, gaint, 'bo', markersize=3, label='auto')
				plt.plot(x, gaints, 'm-', lw=3, label='auto smoothed')
				plt.legend()
				plt.xlim(x.min(), x.max())
				plt.ylim(ymin, ymax)
				plt.xlabel('time', size=16)
				plt.ylabel('[A.U.]', size=16)
			else : 
				plt.plot(x, gaint.real, 'bo', markersize=3, label='real')
				plt.plot(x, gaints.real, 'm-', lw=3, label='real smoothed')
				plt.plot(x, gaint.imag, 'ro', markersize=3, label='imag')
				plt.plot(x, gaints.imag, 'g-', lw=3, label='imag smoothed')
				plt.legend()
				plt.xlim(x.min(), x.max())
				plt.ylim(ymin, ymax)
				plt.xlabel('time', size=16)
				plt.ylabel('[A.U.]', size=16)
			plt.title(r'$G(t) \cdot \sigma_n$, '+strfreq+', '+strbl, size=16)
			plt.savefig(self.outdir+'CaliGain.Gaint_'+strfreq+'_'+strbl.split('=')[0]+'.png')
			if (self.verbose) : print 'CaliGain.Smooth:  end  @', jp.Time(1)+'\n'

