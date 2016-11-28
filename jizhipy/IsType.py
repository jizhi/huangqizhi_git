


class IsType( object ) : 
	dtype = 'class:IsType'

	def isint( self, a ) : 
		typestr = str(type(a)).split("'")[-2].split('.')[-1][:3]
		if (typestr == 'int') : return True
		else : return False
	
	def isfloat( self, a ) : 
		typestr = str(type(a)).split("'")[-2].split('.')[-1][:5]
		if (typestr == 'float') : return True
		else : return False

	def isnum( self, a ) : 
		tf = self.isint(a) + self.isfloat(a)
		if (tf == 0) : return False
		else : return True
	
	def isstr( self, a ) : 
		typestr = str(type(a)).split("'")[-2].split('.')[-1][:3]
		if (typestr == 'str') : return True
		else : return False

	def islist( self, a ) : 
		typestr = str(type(a)).split("'")[-2].split('.')[-1]
		if (typestr == 'list') : return True
		else : return False

	def istuple( self, a ) : 
		typestr = str(type(a)).split("'")[-2].split('.')[-1]
		if (typestr == 'tuple') : return True
		else : return False
	
	def isdict( self, a ) : 
		typestr = str(type(a)).split("'")[-2].split('.')[-1]
		if (typestr == 'dict') : return True
		else : return False
	
	def isndarray( self, a ) : 
		typestr = str(type(a)).split("'")[-2].split('.')[-1]
		if (typestr == 'ndarray') : return True
		else : return False

	def ismatrix( self, a ) : 
		typestr = str(type(a)).split("'")[-2].split('.')[-1]
		if (typestr == 'matrix') : return True
		else : return False

	def ismaskedarray( self, a ) : 
		typestr = str(type(a)).split("'")[-2].split('.')[-1]
		if (typestr == 'MaskedArray') : return True
		else : return False
	
	def isclass( self, a ) : 
		typestr = str(type(a)).split("'")[-2].split('.')[-1]
		if (typestr in ['type', 'classobj']) : return True
		else : return False
		
	def isinstance( self, a ) : 
		typestr1 = str(type(a)).split("'")[-2]
		typestr2 = str(type(a)).split(' ')[0].split('<')[-1]
		if (typestr1=='instance' or typestr2=='class') : 
			return True
		else : return False

