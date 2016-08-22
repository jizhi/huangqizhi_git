from Basic import *
from npfmt import *
from PoolFor import *
import ephem



def _DoMultiprocess_CoordTrans( iterable ) : 
	n1, n2 = iterable[0]
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


def CoordTrans( coordin, coordout, lonORnside, lat=None, epoch='2000', Nprocess=None ) :
	'''
	Transform from coordin ot coordout
	=> 1min transforms 8.3e6 pixels

	coordin, coordout:
		'Galactic' or 'Equatorial'

	lonORnside:
		degree or int(nside)
		longitude or nside
		if (lat is None) : lonORnside = nside of the Healpix Map
		else : lonORnside = longitude (RA or l in Galactic), from 0 to 360 degree

	lat:
		degree
		Dec or b in Galactic, from -90 to +90 degree

	return:
		[lon, lat]
		lon: RA or l, lat: Dec or b in degree
	'''
	#--------------------------------------------------
	coordin, coordout = coordin.lower(), coordout.lower()
	errcoord = False
	if (coordin[0] == 'g') : coordin = 'galactic'
	elif (coordin[0] == 'e') : coordin = 'equatorial'
	else : errcoord = True
	if (coordout[0] == 'g') : coordout = 'galactic'
	elif (coordout[0] == 'e') : coordout = 'equatorial'
	else : errcoord = True
	if (errcoord) : Raise(Exception, 'coordint="'+str(coordout)+'", coordout="'+str(coordout)+'" not in ["galactic","equatorial"]')
	if (coordin == coordout) : return [lonORnside, lat]
	#--------------------------------------------------
	istype, islist = IsType(), False
	if (istype.isint(lonORnside) and lat is None) : 
		n = np.arange(12*lonORnside**2)
		lat, lon = hp.pix2ang(lonORnside, n)
		lat, lon = 90-lat*180/np.pi, lon*180/np.pi
		islist = True
	else : 
		tf = [istype.islist(lonORnside), istype.istuple(lonORnside), istype.isnparray(lonORnside), istype.islist(lat), istype.istuple(lat), istype.isnparray(lat)]
		if (True in tf) : islist = True
		lon = npfmt(lonORnside)
		lat = npfmt(lat)
	if (lon.shape != lat.shape) : Raise(Exception, 'lonORnside.shape='+str(lon.shape)+' != lat.shape='+str(lat.shape))
	#--------------------------------------------------
	shape = lon.shape
	lon, lat = lon.flatten()*np.pi/180, lat.flatten()*np.pi/180
	pool = PoolFor(0, len(lon), Nprocess)
	retn = pool.map_async(_DoMultiprocess_CoordTrans, (lon,lat), (coordin,str(epoch)))
	lon, lat = np.concatenate(retn).T
	lon, lat = lon.reshape(shape), lat.reshape(shape)
	del pool, retn
	#--------------------------------------------------
	if (not islist) : lon, lat = lon.take(0), lat.take(0)
	return [lon*180/np.pi, lat*180/np.pi]



def CoordTransHealpix( healpixmap, coordin, coordout, epoch='2000', Nprocess=None ) : 
	nside = hp.get_nside(healpixmap)
	lon, lat = CoordTrans(coordout, coordin, nside, None, epoch, Nprocess)
	n = hp.ang2pix(nside, (90-lat)*np.pi/180, lon*np.pi/180)
	lon = lat = 0 #2
	return healpixmap[n]

