

class Constant( object ) : 
	dtype = 'class:Constant'
	'''
	Usage:
		const = Constant()
		const.Value('Plank')
		const['Plank']
		const('Plank')
	'''


	def __init__( self ) : 
		self.constant = {'Boltzmann':1.380648813e-23, 
			'Planck':6.6260695729e-34, 
			'c':2.99792458e8, 
			'21cmfreq':1420.40575177e6, 
			'21cmwavelength':0.2110611405413}


	def Value( self, name=None ) : 
		key, value = self.constant.keys(), self.constant.values()
		for i in xrange(len(key)) : key[i] = key[i].lower()
		try : n = key.index(name.lower())
		except : Raise(Exception, 'name='+str(name)+' not in self.constant.keys()='+str(self.constant.keys()))
		return value[n]


	def __call__( self, which='' ) : return self.Value(which)
	def __getitem__( self, which='' ) : return self.Value(which)


