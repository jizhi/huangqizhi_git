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





def _DoMultiprocess_GalacticEquatorial( iterable ) : 
	lon, lat = iterable[1]
	coordin, epoch = iterable[2]
	if (coordin == 'galactic') : 
		func1, func2 = ephem.Galactic, ephem.Equatorial
	else : 
		func1, func2 = ephem.Equatorial, ephem.Galactic
	x, y = [], []
	for i in xrange(len(lon)) : 
		xy = func1(lon[i], lat[i], epoch=epoch)
		xy = func2(xy, epoch=epoch)
		try : x.append(xy.ra+0)
		except : x.append(xy.lon+0)
		try : y.append(xy.dec+0)
		except : y.append(xy.lat+0)
	x, y = npfmt(x), npfmt(y)
	return np.concatenate([x[:,None], y[:,None]], 1)





class CoordTrans( object ) : 


	def __init__( self, Nprocess=None ) : 
		self.Nprocess = Nprocess



	def Gala2Equa( self, l, b, epoch='2000' ) : 
		'''
		l, b: 
			in degree
			l, b can be one or N-D array (any shape)
			l and b can broadcast !
		'''
		l, b, islist = self._CheckList(l, b)
		return self._GalacticEquatorial(l, b, 'galactic', epoch, islist)



	def Equa2Gala( self, RA, Dec, epoch='2000' ) : 
		'''
		RA, Dec:
			in degree
			RA, Dec can be one or N-D array (any shape)
			RA and Dec can broadcast !
		'''
		RA, Dec, islist = self._CheckList(RA, Dec)
		return self._GalacticEquatorial(RA, Dec, 'equatorial', epoch, islist)



	def Gala2EquaHeal( self, healpixmap, ordering='RING', epoch='2000' ):
		'''
		healpixmap:
			Healpix map with nside=2**n, total pixels=12*nside**2

		ordering: 
			'RING' or 'NESTED'
		'''
		nside = hp.get_nside(healpixmap)
		RA, Dec, nest = self._Nside2LonLat(nside, ordering)
		l, b = self.Equa2Gala(RA, Dec, epoch)
		n = hp.ang2pix(nside, (90-b)*np.pi/180, l*np.pi/180, nest=nest)
		return healpixmap[n]



	def Equa2GalaHeal( self, healpixmap, ordering='RING', epoch='2000' ):
		'''
		healpixmap:
			Healpix map with nside=2**n, total pixels=12*nside**2

		ordering: 
			'RING' or 'NESTED'
		'''
		nside = hp.get_nside(healpixmap)
		l, b, nest = self._Nside2LonLat(nside, ordering)
		RA, Dec = self.Gala2Equa(l, b, epoch)
		n =hp.ang2pix(nside, (90-Dec)*np.pi/180, RA*np.pi/180, nest=nest)
		return healpixmap[n]



	def Celestial( self, RAl, Decb, coordin, coordout, epoch='2000' ) : 
		'''
		RAl, Decb:
			in degree
			Can be one or N-D array
		'''
		coordin, coordout = self._CoordInOut(coordin, coordout)
		if (coordin == coordout) : return np.array([RAl, Decb])
		elif (coordin == 'galactic') : return self.Gala2Equa(RAl, Decb, epoch)
		elif (coordin == 'equatorial') : return self.Equa2Gala(RAl, Decb, epoch)

		

	def CelestialHeal( self, healpixmap, coordin, coordout, ordering='RING', epoch='2000' ) : 
		coordin, coordout = self._CoordInOut(coordin, coordout)
		if (coordin == coordout) : return healpixmap
		elif (coordin == 'galactic') : return self.Gala2EquaHeal(healpixmap, ordering, epoch)
		elif (coordin == 'equatorial') : return self.Equa2GalaHeal(healpixmap, ordering, epoch)





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
			in degree
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







	def _CheckList( self, x, y ) : 
		istype = IsType()
		if (istype.isnum(x) and istype.isnum(y)) : islist = False
		else : islist = True
		x, y = npfmt(x), npfmt(y)
		try : 
			x = x + y*0
			y = y + x*0
		except : Raise(Exception, 'x.shape='+str(x.shape)+', y.shape='+str(y.shape)+', can NOT broadcast')
		return [x, y, islist]





	def _GalacticEquatorial( self, lon, lat, coordin, epoch, islist ) : 
		shape = lon.shape
		lon, lat = lon.flatten()*np.pi/180, lat.flatten()*np.pi/180
		pool = PoolFor(0, len(lon), self.Nprocess)
		retn = pool.map_async(_DoMultiprocess_GalacticEquatorial, (lon,lat), (coordin, str(epoch)))
		lon, lat = np.concatenate(retn).T
		lon, lat = lon.reshape(shape), lat.reshape(shape)
		if (not islist) : lon, lat = lon.take(0), lat.take(0)
		return np.array([lon*180/np.pi, lat*180/np.pi])  # degree



	def _Nside2LonLat( self, nside, ordering ) : 
		''' return in degree '''
		ordering = str(ordering).lower()
		if (ordering == 'ring') : nest = False
		else : nest = True
		lat, lon = hp.pix2ang(nside, np.arange(12*nside**2), nest=nest)
		lat = np.pi/2 - lat
		return [lon*180/np.pi, lat*180/np.pi, nest]



	def _Angle( self, kwargs, orderkwargs ) : 
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





	def Heal2Flat( self, lon, lat, thetawidth, Npix, hpmapORnside, ordering='RING' ) : 
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



