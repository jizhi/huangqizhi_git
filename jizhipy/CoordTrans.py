from Other import *
from npfmt import *
from PoolFor import *
from IsType import *
from OrderKwargs import *
from ArrayAxis import *
import ephem
import healpy as hp
from Print import *
from BeamModel import *




def _DoMultiprocess_Celestial( iterable ) : 
	RAl, Decb = iterable[1]
	coordin, coordout, epochin, epochout = iterable[2]
	if   (coordin  == 'galactic')   : funcin  = ephem.Galactic
	elif (coordin  == 'equatorial') : funcin  = ephem.Equatorial
	if   (coordout == 'galactic')   : funcout = ephem.Galactic
	elif (coordout == 'equatorial') : funcout = ephem.Equatorial
	x, y = [], []
	for i in xrange(len(RAl)) : 
		xy = funcin(RAl[i], Decb[i], epoch=epochin)
		xy = funcout(xy, epoch=epochout)
		if   (coordout == 'galactic') : 
			x.append(xy.lon +0)
			y.append(xy.lat +0)
		elif (coordout == 'equatorial') : 
			x.append(xy.ra  +0)
			y.append(xy.dec +0)
	return np.array([x, y])





class CoordTrans( object ) : 


	def __init__( self, Nprocess=None ) : 
		self.Nprocess = Nprocess


	######## (coordin, epochin) => (coordout, epochout) ########


	def _CoordInOut( self, coordin, coordout ) : 
		raisestr = "coordin="+str(coordin)+", coordout="+str(coordout)+" not in ['galactic', 'equatorial']"
		coordin, coordout = str(coordin).lower(), str(coordout).lower()
		if   ('gal' == coordin[:3]) : coordin = 'galactic'
		elif ('equ' == coordin[:3]) : coordin = 'equatorial'
		else : Raise(Exception, raisestr)
		if   ('gal' == coordout[:3]) : coordout = 'galactic'
		elif ('equ' == coordout[:3]) : coordout = 'equatorial'
		else : Raise(Exception, raisestr)
		return [coordin, coordout]



	def _EpochInOut( self, epochin, epochout ) : 
		'''
		epoch:
			(1) '2000' or ephem.J2000
			(2) '1950' or ephem.B1950
			(3) '1900' or ephem.B1900
			(4) other number
		'''
		def _Epoch( epoch ) : 
			istype = IsType()
			if (istype.isstr(epoch) or istype.isnum(epoch)) : 
				epoch = str(epoch)
				if   ('2000' in epoch) : epoch = ephem.J2000
				elif ('1950' in epoch) : epoch = ephem.B1950
				elif ('1900' in epoch) : epoch = ephem.B1900
				else : epoch = str(epoch)
			return epoch
		return [_Epoch(epochin), _Epoch(epochout)]



	def _RAlDecb( self, RAl, Decb ) : 
		istype = IsType()
		if (istype.isnum(RAl) and istype.isnum(Decb)) : 
			islist = False
		else : islist = True
		RAl, Decb = npfmt(RAl), npfmt(Decb)
		try : 
			RAl  = RAl + Decb*0
			Decb = RAl*0 + Decb
		except : Raise(Exception, 'RAl.shape='+str(RAl.shape)+', Decb.shape='+str(Decb.shape)+', can NOT broadcast')
		return [RAl, Decb, islist]




	def Celestial( self, RAl, Decb, coordin, coordout, epochin='2000', epochout='2000' ) : 
		'''
		RAl, Decb:
			in rad, NOT degree
			Can be any shape
			But must can broadcast

		return:
			[RAl, Decb]
				in rad
				with same shape as (RAl+Decb).shape
		'''
		coordin, coordout = self._CoordInOut(coordin, coordout)
		epochin, epochout = self._EpochInOut(epochin, epochout)
		if (coordin==coordout and epochin==epochout) : 
			return [RAl, Decb]
		RAl, Decb, islist = self._RAlDecb(RAl, Decb)
		shape = RAl.shape
		sent = [RAl.flatten(), Decb.flatten()]
		bcast = [coordin, coordout, epochin, epochout]
		if (self.Nprocess <= 1) : 
			x, y = _DoMultiprocess_Celestial([None, sent, bcast])
		else : 
			pool = PoolFor(0, len(RAl), self.Nprocess)
			xy=pool.map_async(_DoMultiprocess_Celestial, sent,bcast)
			x, y = np.concatenate(xy, 1)
		x, y = x.reshape(shape), y.reshape(shape)
		if (not islist) : x, y = x[0], y[0]
		return np.array([x, y])


	#### Healpix, (coordin, epochin) => (coordout, epochout) ####


	def _Nside( self, nside ) : 
		n = np.log(nside) / np.log(2)
		if (n != int(n)) : Raise(Exception, 'nside='+str(nside)+' != 2**n')
		return int(round(nside))



	def _Ordering( self, ordering ) : 
		ordering = str(ordering).upper()
		if ('NEST' in ordering) : ordering, nest = 'NESTED', True
		else : ordering, nest = 'RING', False
		return [ordering, nest]




	def CelestialHealpix( self, healpixmap, ordering, coordin, coordout, epochin='2000', epochout='2000' ) : 
		'''
		healpixmap:
			(1) np.ndarray with shape=(12*nside**2,)
			(2) int number =nside

		ordering:
			'RINGE' | 'NESTED'

		return:
			Case (1) [hpix_in2out, healpixmap_in2out]
			Case (2) hpix_in2out

				Usage: healpixmap_in2out = healpixmap[hpix_in2out]
		'''
		try : nside = hp.get_nside(healpixmap)
		except : nside = self._Nside(healpixmap)
		ordering, nest = self._Ordering(ordering)
		hpix_in2out = np.arange(12*nside**2)
		Decb, RAl = hp.pix2ang(nside, hpix_in2out, nest=nest)
		Decb = np.pi/2 - Decb
		RAl, Decb = self.Celestial(RAl, Decb, coordout, coordin, epochout, epochin)
		hpix_in2out =hp.ang2pix(nside, np.pi/2-Decb,RAl, nest=nest) 
		try : return [hpix_in2out, healpixmap[hpix_in2out]]
		except : return hpix_in2out


	################# xyz, theta, Rotation #################


	def _Angle( self, kwargs, orderkwargs ) : 
		''' degree to rad '''
		N = len(orderkwargs)
		islist, istype = [], IsType()
		raisestr = 'Shape miss-matching, '
		for i in xrange(N) : 
			ok = orderkwargs[i]
			islist.append( 1-istype.isnum(kwargs[ok]) )
			kwargs[ok] = npfmt(kwargs[ok]) *np.pi/180
			raisestr += ok+'.shape='+str(kwargs[ok].shape)+', '
		try : 
			for i in xrange(N) : 
				oki = orderkwargs[i]
				for j in xrange(N) : 
					okj = orderkwargs[j]
					if (i == j) : continue
					kwargs[oki] = kwargs[oki] + kwargs[okj]*0  # broadcast
		except : jp.Raise(Exception, raisestr)
		if (np.array(islist).sum() == 0) : islist = False
		else : islist = True
		return [kwargs, islist]



	def _RotationMatrix( self, key, ang, islist ) : 
		one  = np.ones(ang.shape)
		zero = np.zeros(ang.shape)
		if (key == 'ax') : 
			R = np.array([[ one,     zero    ,    zero    ],
			              [zero,  np.cos(ang), np.sin(ang)],
			              [zero, -np.sin(ang), np.cos(ang)]])
		elif (key == 'ay') : 
			R = np.array([[np.cos(ang), zero, -np.sin(ang)],
			              [   zero   ,   one,     zero    ],
			              [np.sin(ang), zero,  np.cos(ang)]])
		elif (key == 'az') : 
			R = np.array([[ np.cos(ang), np.sin(ang), zero],
			              [-np.sin(ang), np.cos(ang), zero],
			              [    zero    ,    zero    ,  one]])
		else : Raise(Exception, "key in **kwargs not in ['ax', 'ay', 'az']")
		R = ArrayAxis(R, 0, -1, 'move')
		R = ArrayAxis(R, 0, -1, 'move')
		if (not islist) : R = R[0]
		return R




	def xyzRotationMatrix( self, **kwargs ) : 
		'''
		Note that here rotate the xyz coordinate and the point doesn't move, NOT rotate the point

		Right-hand rule:
			Your right thumb points along the +Z axis and the curl of your fingers represents a motion from +X to +Y to -X to -Y. When viewed from the top along -Z, the system is counter-clockwise.
			The angle along the fingers is positive.

		**kwargs:
			keywords must in [ax, ay, az, order]
			Order of the keywords is the order of the rotations
		
		xyzRotationMatrix(ay, ax, az) : 
			first rotate around Y-axis, second rotate around X-axis, third rotate around Z-axis

		ax: Rotation angle around X-axis, thumb points to +X
		ay: Rotation angle around Y-axis, thumb points to +Y
		az: Rotation angle around Z-axis, thumb points to +Z
		order: give the order manually, ['az', 'ax', 'ay']

		ax, ay, az:
			in degree, NOT rad
			Can be one or N-D array

		return:
			Rotation matrix R in np.array(), NOT np.matrix()
				new_xyz = R * old_xyz

			R.shape = ax.shape+(3,3)
		'''
		try : orderkwargs = kwargs['order']
		except : orderkwargs = OrderKwargs(2)
		N = len(orderkwargs)
		kwargs, islist = self._Angle(kwargs, orderkwargs)
		R = []
		for i in xrange(N) : 
			ok = orderkwargs[i]
			Ri = self._RotationMatrix(ok, kwargs[ok], islist)
			if (not islist) : Ri = Ri[None,:]
			shape = Ri.shape  #@
			Ri = Ri.reshape(np.prod(shape[:-2]), 3, 3)  # (N,3,3)
			R.append(Ri)
		for i in xrange(1, N) : 
			for j in xrange(len(R[0])) : 
				R[0][j] = np.array(np.matrix(R[i][j]) * np.matrix(R[0][j]))
		R = R[0]
		R = R.reshape(shape)
		if (not islist) : R = R[0]
		return R





	def xyzRotation( self, xyz, **kwargs ) : 
		'''
		xyz:
			xyz can be any shape, but must:
			x, y, z = xyz

			x = sin(theta) * cos(phi)
			y = sin(theta) * sin(phi)
			z = cos(theta)

		**kwargs:
			See self.xyzRotationMatrix()

		return:
			Same shape and type as input xyz

			xyz_new.shape = xyz.shape+(3,)   # (3,) for x,y,z
		'''
		xyz = np.array(xyz, float).T
		if (xyz.shape == (3,)) : 
			islistx = False
			xyz = xyz[None,:]
		else : islistx = True
		shapex = xyz.shape
		if (shapex[-1] != 3) : Raise(Exception, 'xyz.shape='+str(shapex[::-1])+', shape[0] != 3')
		xyz = xyz.reshape(np.prod(shapex[:-1]), 3)
		x, y, z = xyz.T
		#--------------------------------------------------
		try : kwargs['order']
		except : kwargs['order'] = OrderKwargs()
		R = self.xyzRotationMatrix(**kwargs)
		if (R.shape == (3,3)) : 
			R = R[None,:]
			islistr = False
		else : islistr = True
		shaper = R.shape
		R = R.reshape(np.prod(shaper[:-2]), 3, 3)
		Rx, Ry, Rz = R[:,:,0], R[:,:,1], R[:,:,2]  # (M,3)
		#--------------------------------------------------
		xyz = x[:,None,None]*Rx[None,:] + y[:,None,None]*Ry[None,:] + z[:,None,None]*Rz[None,:]
		shapex, shaper = list(shapex[:-1]), list(shaper[:-2])
		if (not islistx) : shapex = []
		if (not islistr) : shaper = []
		xyz = xyz.reshape(shapex + shaper + [3])
		return xyz.T





	def thetaphi2xyz( self, thetaphi ) : 
		'''
		thetaphi:
			in rad, NOT degree
			Can be any shape, but must:
			theta, phi = thetaphi

			x = sin(theta) * cos(phi)
			y = sin(theta) * sin(phi)
			z = cos(theta)

		return:
			xyz = [x, y, z]
			x.shape = y.shape = z.shape = theta.shape = phi.shape
		'''
		theta, phi = thetaphi
		if (IsType().isnum(theta)) : islisttheta = False
		else : islisttheta = True
		if (IsType().isnum(phi)) : islistphi = False
		else : islistphi = True
		islist = bool(islisttheta + islistphi)
		theta, phi = npfmt(theta), npfmt(phi)
		if (theta.shape != phi.shape) : Raise(Exception, 'theta.shape='+str(theta.shape)+' != phi.shape='+str(phi.shape))
		x = np.sin(theta) * np.cos(phi)
		y = np.sin(theta) * np.sin(phi)
		z = np.cos(theta)
		xyz = np.array([x, y, z])
		x = y = z = 0 #@
		if (not islist) : xyz = xyz.flatten()
		return xyz





	def xyz2thetaphi( self, xyz ) : 
		'''
		xyz:
			Can be any shape, but must:
			x, y, z = xyz

		return:
			thetaphi = [theta, phi]
			in rad (NOT degree)
			theta.shape = phi.shape = x.shape
		'''
		x, y, z = xyz
		if (IsType().isnum(x)) : islistx = False
		else : islistx = True
		if (IsType().isnum(y)) : islisty = False
		else : islisty = True
		if (IsType().isnum(z)) : islistz = False
		else : islistz = True
		islist = bool(islistx + islisty + islistz)
		x, y, z = npfmt(x), npfmt(y), npfmt(z)
		#--------------------------------------------------
		# theta from 0 to 180, sin(theta) >=0
		theta = np.arccos(z)  # 0 to 180 degree
		zero = (np.sin(theta) == 0)
		sign = np.sign(z[zero])
		theta[zero] += 1e-30
		del z
		WarningFilter(True)
		x = x / np.sin(theta)
		y = y / np.sin(theta)
		WarningFilter(None)
		x = x + 1j*y
		del y
		phi = np.angle(x) % (2*np.pi)
		theta[zero] -= 1e-30
		if (zero.sum() > 0) : phi[zero] = np.pi/2 - sign*np.pi/2	
		return np.array([theta, phi])





	def thetaphiRotation( self, thetaphi, **kwargs ) : 
		'''
		thetaphi:
			in rad, NOT degree
			thetaphi can be any shape, but must:
			theta, phi = thetaphi
		'''
		theta, phi = np.array(thetaphi, float)
		shape = theta.shape
		xyz = np.array([np.sin(theta)*np.cos(phi), np.sin(theta)*np.sin(phi), np.cos(theta)])
		del theta, phi
		kwargs['order'] = OrderKwargs()
		xyz = self.xyzRotation(xyz, **kwargs)
		thetaphi = self.xyz2thetaphi(xyz)
		return thetaphi


	##################################################


	def Healpix2Flat( self, lon, lat, thetawidth, Npix, hpmapORnside, ordering='RING' ) : 
		'''
		lon, lat:
			in degree
			lon=RA/l, lat=Dec/b
			Center of the region
	
		thetawidth:
			in degree
			Width of the region (theta angle)

		hpmapORnside:
			(1) is nside: return [thetaphi, healpix_index]
			(2) is healpixmap: return [thetaphi, healpix_index, sky_region]

			NOTE THAT thetaphi[0] is theta = 90-lat

		return:
			[ (theta,phi), hppix ]
		 OR
			[ (theta,phi), hppix, hpmapORnside[hppix] ]
		'''
		lon, lat, thetawidth = Num(lon, float)%360, Num(lat, float), Num(thetawidth, float)*np.pi/180
		if (str(ordering).lower() == 'ring') : nest = False
		else : nest = True
		if (IsType().isnum(hpmapORnside)) : nside = hpmapORnside
		else : nside = hp.get_nside(hpmapORnside)
		theta, phi = ThetaPhiMatrix(thetawidth, Npix)  # rad
		# Around x-axis, rotate lon
		# Then around y-axis, rotate -(90-lat)
		theta, phi = self.thetaphiRotation([theta, phi], ay=-(90-lat), az=-lon)
		theta, phi = theta.T, phi.T[::-1,::-1]
		hppix = hp.ang2pix(nside, theta, phi, nest=nest)
		if (phi.max()-phi.min() > 2*np.pi-thetawidth) : 
			phi = np.pi - phi
			if ((abs(phi)<2*thetawidth).sum() != 0) : 
				phi = np.pi - phi
			else : 
				phi[phi<0] = abs(phi[phi<0]) - np.pi
				phi[phi>=0] = np.pi - phi[phi>=0]
		if (IsType().isnum(hpmapORnside)) : return [np.array([theta,phi]), hppix]
		else : return [np.array([theta,phi]), hppix, hpmapORnside[hppix]]



