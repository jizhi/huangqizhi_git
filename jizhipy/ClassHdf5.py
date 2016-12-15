from IsType import *
from Path import *
import numpy as np
import h5py



class ClassHdf5( object ) : 


	def __init__( self, classinstance, outname, verbose=True ) : 
		'''
		classinstance:
			Instance of class: a = A()

		outname:
			end with '.hdf5'
		'''
		self.classinstance = classinstance
		self.outname = outname
		self.verbose = bool(verbose)
		self.keys, self.value = None, None



	def Keys( self ) : 
		'''
		return:
			all keys in a
		'''
		keys = self.classinstance.__dict__.keys()
		istype = IsType()
		dowhile, times = True, 0
		while(dowhile) : 
			dowhile, times = False, times+1
			for i in xrange(len(keys)) : 
				k = keys[i].split('.')
				v = self.classinstance
				for j in xrange(len(k)) : v = v.__dict__[k[j]]
				if (istype.isfunc(v) or istype.isclass(v)) : 
					keys[i] = ''
				elif (istype.isinstance(v) and not istype.isspecialinstance(v)) : 
					k = v.__dict__.keys()
					for j in xrange(len(k)) : k[j] = keys[i] + '.'+k[j]
					keys[i] = k
					dowhile = True
			keytmp = []
			for i in xrange(len(keys)) : 
				if (keys[i] == '') : continue
				elif (type(keys[i]) == list) : keytmp += keys[i]
				else : keytmp.append(keys[i])
			keys = keytmp
		self.keys = keys





	def Values( self ) : 
		'''
		return:
			All values corresponds to keys
		'''
		if (self.keys is None) : self.Keys()
		self.values = []
		for i in xrange(len(self.keys)) : 
			k = self.keys[i].split('.')
			v = self.classinstance
			for j in xrange(len(k)) : 
				v = v.__dict__[k[j]]
			self.values.append( v )





	def Save( self ) : 
		'''
		Save instance to .hdf5
		'''
		self.outname = str(self.outname)
		if (self.outname[-5:] !='.hdf5') : self.outname += '.hdf5'
		ExistsPath(self.outname, old=True)
		if (self.verbose) : print 'jizhipy.ClassHdf5: Saveing to  "'+self.outname+'"'
		fo = h5py.File(self.outname, 'w')
		if (self.value is None or self.keys is None): self.Values()
		istype = IsType()
		fo['classinstance'] = str(type(self.classinstance))
		for i in xrange(len(self.keys)) : 
			if (self.values[i] is None) : 
				fo[self.keys[i]] = 'None'
				fo[self.keys[i]].attrs['type'] = 'is None, not str'
			#--------------------------------------------------
			elif (istype.ismatrix(self.values[i])) : 
				fo[self.keys[i]] = np.array(self.values[i])
				fo[self.keys[i]].attrs['type'] = 'numpy.matrix'
			#--------------------------------------------------
			elif (istype.ismaskedarray(self.values[i])) : 
				fo[self.keys[i]] = self.values[i].data
				fo[self.keys[i]+'.mask'] = self.values[i].mask
				fo[self.keys[i]].attrs['type'] = 'numpy.MaskedArray .data'
				fo[self.keys[i]+'.mask'].attrs['type'] = 'numpy.MaskedArray .mask'
			#--------------------------------------------------
			elif (istype.isdict(self.values[i])) : 
				fo[self.keys[i]] = self.values[i].values()
				fo[self.keys[i]+'.keys'] = self.values[i].keys()
				fo[self.keys[i]].attrs['type'] = 'dict .values()'
				fo[self.keys[i]+'.keys'].attrs['type'] = 'dict .keys()'
			#--------------------------------------------------
			else : 
				fo[self.keys[i]] = self.values[i]
				if   (istype.islist(self.values[i])) : 
					fo[self.keys[i]].attrs['type'] = 'list'
				elif (istype.istuple(self.values[i])) : 
					fo[self.keys[i]].attrs['type'] = 'tuple'
				else : fo[self.keys[i]].attrs['type'] = 'numpy.array'
			#--------------------------------------------------
		fo.flush()
		fo.close()





	def Read( self ) : 
		'''
		When read .hdf5 to a class:instance, it needs this class:instance already has all necessary sub-instance !
		'''
		self.outname = str(self.outname)
		if (self.outname[-5:] !='.hdf5') : self.outname += '.hdf5'
		if (self.verbose) : print 'jizhipy.ClassHdf5: Reading from  "'+self.outname+'"'
		fo = h5py.File(self.outname, 'r')
		keys = fo.keys()
		for i in xrange(len(keys)) : 
			keys[i] = str(keys[i])
			if (keys[i]=='classinstance' or keys[i][-5:]=='.mask' or keys[i][-5:]=='.keys') : continue
			Type = fo[keys[i]].attrs['type']
			k = keys[i].split('.')
			v = self.classinstance
			for j in xrange(len(k)-1) : v = v.__dict__[k[j]]
			#--------------------------------------------------
			if ('None' in Type) : 
				v.__dict__[k[-1]] = None
			#--------------------------------------------------
			elif ('numpy.matrix' in Type) : 
				v.__dict__[k[-1]] = np.matrix(fo[keys[i]].value)
			#--------------------------------------------------
			elif ('numpy.MaskedArray' in Type) : 
				value = np.ma.MaskedArray(fo[keys[i]].value)
				value.mask = fo[keys[i]+'.mask'].value
				v.__dict__[k[-1]] = value 
			#--------------------------------------------------
			elif ('dict' in Type) : 
				fokey = fo[keys[i]+'.keys'].value
				fovalue = fo[keys[i]].value
				value = {}
				for j in xrange(len(fokey)) : 
					value[fokey[j]] = fovalue[j]
				v.__dict__[k[-1]] = value
			#--------------------------------------------------
			elif (Type == 'list') : 
				v.__dict__[k[-1]] = list(fo[keys[i]].value)
			#--------------------------------------------------
			elif (Type == 'tuple') : 
				v.__dict__[k[-1]] = tuple(fo[keys[i]].value)
			#--------------------------------------------------
			else : v.__dict__[k[-1]] = fo[keys[i]].value
		#--------------------------------------------------
		fo.close()
		return self.classinstance









#
#class A(object) : 
#	def __init__(self) : 
#		self.a1, self.a2 = None, (1,2)
#
#class B(object) : 
#	def __init__(self) : 
#		self.b1, self.b2, self.a = ['b1','b2'], np.arange(2), A()
#
#class C(object) : 
#	def __init__(self) : 
#		self.c1, self.c2, self.a, self.b = np.matrix([[1,2],[3,4]]), np.ma.MaskedArray(np.arange(2), [True,False]), A(), B()
#
#class D(object) : 
#	def __init__(self) : 
#		self.d1, self.d2, self.a, self.b, self.c = 'whoami', 789, A(), B(), C()
#
#
#
#outname = 'test.hdf5'
#
#a = D()
#classhdf5 = jp.ClassHdf5(a, outname)
#classhdf5.Save()
#print a.a.a1, a.a.a2
#print a.b.b1, a.b.b2, a.a.a1, a.a.a2
#print a.c.c1, a.c.c2, a.c.b.b1, a.c.b.b2, a.c.b.a.a1, a.c.b.a.a2
#print a.d1, a.d2
#print
#print '---------------'
#print
#
#
#b = D()
#b.a.a1, b.a.a2 = None, None
#b.b.b1, b.b.b2, b.a.a1, b.a.a2 = None, None, None, None
#b.c.c1, b.c.c2, b.c.b.b1, b.c.b.b2, b.c.b.a.a1, b.c.b.a.a2 = None, None, None, None, None, None
#b.d1, b.d2 = None, None
#
#classhdf5 = jp.ClassHdf5(b, outname)
#classhdf5.Read()
#
#print b.a.a1, b.a.a2
#print b.b.b1, b.b.b2, b.a.a1, b.a.a2
#print b.c.c1, b.c.c2, b.c.b.b1, b.c.b.b2, b.c.b.a.a1, b.c.b.a.a2
#print b.d1, b.d2
#


