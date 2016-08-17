from CaliPhase import *
from CaliTemperature import *
##################################################


class MapMaking( object ) : 


	def __init__( self, caliphase=None, calitemp=None ) : 
		if (caliphase is not None) : self.caliphase = caliphase
		if (calitemp is not None) : 
			self.calitemp = calitemp
			self.vis = self.calitemp.vis
			self.calitemp.__dict__.pop('vis') # will also pop('vis') of the calitemp outside
			self._Outdir()


	def _Outdir( self ) : 
		outdir = sys.argv[0][:-3] + '_output/'
		if (outdir[:2]=='./') : outdir = outdir[2:]
		outdir += self.calitemp.antarray.Hdf5.hdf5dir[:-1].split('/')[-1]+'/'
		self.outdir = outdir + self.dtype[6:]+'/'
		mkdir(self.outdir)
		

	def CaliPhase( self ) : 
		nfreq = self.calitemp.nfreq
		try : Padd = self.caliphase.Phaseaddv[nfreq]
		except AttributeError : Padd = self.caliphase.Phaseaddp[nfreq]
		except : Padd = self.caliphase.Phaseadde[nfreq]
		Pns = self.caliphase.Phasens[nfreq,:][None,:]
		Padd = Padd[None,:]
		nauto = self.calitemp.antarray.Blorder.auto.size
		self.vis.data[:,nauto:] *= np.exp(-1j*(Pns + Padd))


	def MapMaking( self, Nangle=200 ) : 
		''' anglesize: in rad '''
		nauto = self.calitemp.antarray.Blorder.auto.size
		mask = self.vis.mask.sum(1)
		auto = self.vis.data[:,:nauto].sum(1).real 
		cross = self.vis.data[:,nauto:].sum(1).real
		self.mapmaking = auto + 2*cross
		self.mapmaking /= nauto  # real ndarray
		np.save(self.outdir+'mapmaking.data.npy', self.mapmaking)
		np.save(self.outdir+'mapmaking.mask.npy', mask, bool)
		#--------------------------------------------------
		visorder = self.caliphase.antarray.visorder
		bl = self.caliphase.antarray.Blorder.blorder[visorder]
		blmax = abs(bl[:,0]).max()
		fwhm = 1.03*300/self.calitemp.freq/blmax
		angle = self._MapMaking1to2(fwhm, Nangle)
		mask = np.zeros(self.mapmaking.shape,bool) +mask[None,:]
		self.mapmaking = np.ma.MaskedArray(self.mapmaking, mask)  # (200,Nt)
		mask = 0 #@
		#--------------------------------------------------
		# timem to RA
		shape = self.caliphase.antarray.vis.shape
		bl = abs(bl[:,0])
		nv = np.where(bl==bl.max())[0][0]
		nv = visorder[nv]
		vis = abs(self.caliphase.antarray.vis[:,shape[1]/2,nv])
		nt = np.where(vis==vis.max())[0][0]
		tmin = (self.caliphase.antarray.Hdf5.nhdf5*shape[0] + nt) * self.caliphase.antarray.Ant.inttime /60.
		RA = self.caliphase.RA
		self.RA = (RA + (self.calitemp.timem - tmin)/60.*15) %360
		if (358 < self.RA[-1]-self.RA[0] <= 360) : 
			for n in xrange(1, len(self.RA)) : 
				if (self.RA[n]-self.RA[n-1] < 350) : break
			self.RA = np.concatenate([self.RA[n:], self.RA[:n]])
			self.mapmaking = np.concatenate(self.mapmaking[:,n:], self.mapmaking[:,:n], 1)
		self.Dec = self.caliphase.Dec + angle *180/np.pi
		np.save(self.outdir+'RAlist.npy', self.RA)
		np.save(self.outdir+'Declist.npy', self.Dec)
		#--------------------------------------------------
		plt.plot(self.RA, self.mapmaking.data[Nangle/2])
		plt.xlim(self.RA.min(), self.RA.max())
		plt.xlabel('RA [deg]', size=16)
		plt.ylabel('T [K]', size=16)
		plt.title('1D synthesic map @ '+str(self.calitemp.freq)+'MHz', size=16)
		plt.savefig(self.outdir+'mapmaking_1D.png')
		plt.close()
		#--------------------------------------------------
		plt.figure(figsize=(9,4))
		plt.pcolormesh(self.RA, self.Dec, self.mapmaking.data)
		cbar = plt.colorbar()
		cbar.set_label('K', fontsize=16)
		plt.xlim(self.RA.min(), self.RA.max())
		plt.xlabel('RA [deg]', size=16)
		plt.ylim(self.Dec.min(), self.Dec.max())
		plt.ylabel('Dec [deg]', size=16)
		plt.title('Maked map @ '+str(self.calitemp.freq)+'MHz', size=16)
		plt.savefig(self.outdir+'mapmaking_2D.png')
		plt.close()
		

	def _MapMaking1to2( self, fwhm, Nangle ) : 
		if (Nangle <= 1) : return
		anglesize = 2*fwhm
		N = Nangle
		if (N%2 == 0) : N +=1
		angle = anglesize/fwhm * 2*np.pi/3
		angle = np.linspace(-angle/2, angle/2, N)
		tf = abs(angle)>np.pi/2
		angle[:N/2][tf[:N/2]] = -np.pi/2
		angle[N/2:][tf[N/2:]] = np.pi/2
		angle = angle[:Nangle]
		self.mapmaking = self.mapmaking[None,:] * np.cos(angle[:,None])  # (200,Nt)
		np.save(self.outdir+'mapmaking1to2.npy', np.cos(angle))
		angle = np.linspace(-fwhm, fwhm, N)[:Nangle]
		return angle


