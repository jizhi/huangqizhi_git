from Raise import *
from npfmt import *



def Convolve( source, beam, edge='mirror' ) : 
	'''
	This function is used to convolve two nd-array with same shape (source.shape==beam.shape).

	beam:
		You can normalize it to beam.sum()=1 or beam.max()=1 or not normalize depending on the situation.
	
	edge: 
		How to handle the edge.
		(1) edge='linear', convolve directly.
		(2) edge=value(0, 1, source.min() and so on), outside of the source will be treated as this value.
		(3) edge='mirror', the edge as a mirror.
	'''
	source, beam = npfmt(source), npfmt(beam)
	shape, dim = source.shape, len(source.shape)
	edge = str(edge).lower()
	if (source.shape != beam.shape) : 
		Raise(Exception, 'source.shape='+str(source.shape)+' != beam.shape='+str(beam.shape))
	if (dim > 3) : Raise(Exception, 'source/beam is '+str(dim)+'D, but now can just handle <=3D. You can modify this function yourself for >=4D')

	def fft( x ) : 
		if   (dim == 1) : y = np.fft.fft( x ) 
		elif (dim == 2) : y = np.fft.fft2( x ) 
		else : y = np.fft.fftn( x ) 
		return y

	def ifft( x ) : 
		if   (dim == 1) : y = np.fft.ifft( x ) 
		elif (dim == 2) : y = np.fft.ifft2( x ) 
		else : y = np.fft.ifftn( x ) 
		return y

	# npix should be even when convolve (FFT)
	# Modify here for >=4D
	if (shape[0]%2 == 1) : 
		source = np.append(source, source[-1:], 0)
		beam   = np.append(  beam,   beam[-1:], 0)
	if (dim >= 2) : 
		if (shape[1]%2 == 1) : 
			source = np.append(source, source[:,-1:], 1)
			beam   = np.append(  beam,   beam[:,-1:], 1)
	if (dim >= 3) : 
		if (shape[2]%2 == 1) : 
			source = np.append(source, source[:,:,-1:], 2)
			beam   = np.append(  beam,   beam[:,:,-1:], 2)

	if (edge == 'linear') : 
		image = np.fft.fftshift(ifft(fft(source)*fft(beam))).real
	else :
		if (edge == 'mirror') : v, edge = 1, 0
		else : v, edge = 0, float(edge)
		
		if (dim == 1) : 
			n0 = source.shape[0]
			n1 = int(round(n0/2.))
			source1 = source[::-1] *v+edge
			source = np.append(source1, source)
			source = np.append(source, source1)
			source = source[n0-n1:n1-n0]
			source1 = 0 #@
			# For beam, always add 0
			beam1 = np.zeros([n0+2*n1,], beam.dtype)
			beam1[n0-n1:n1-n0] = beam
			beam = beam1
			image = np.fft.fftshift(ifft(fft(source)*fft(beam))).real
			image = image[n1:-n1]

		elif (dim == 2) : 
			n0 = npfmt(source.shape)
			n1 = (n0/2.).round().astype(int)
			source1 = source[:,::-1] *v+edge
			source = np.append(source1, source, 1)
			source = np.append(source, source1, 1)
			source1 = source[::-1]
			source = np.append(source1, source, 0)
			source = np.append(source, source1, 0)
			source = source[n0[0]-n1[0]:n1[0]-n0[0],n0[1]-n1[1]:n1[1]-n0[1]]
			source1 = 0 #@
			# For beam, always add 0
			beam1 = np.zeros(n0+2*n1, beam.dtype)
			beam1[n0[0]-n1[0]:n1[0]-n0[0],n0[1]-n1[1]:n1[1]-n0[1]] = beam
			beam = beam1
			image = np.fft.fftshift(ifft(fft(source)*fft(beam))).real
			image = image[n1[0]:-n1[0],n1[1]:-n1[1]]

		elif (dim == 3) : 
			n0 = npfmt(source.shape)
			n1 = (n0/2.).round().astype(int)
			# axis 0
			source1 = source[::-1] *v+edge
			source = np.append(source1, source, 0)
			source = np.append(source, source1, 0)
			# axis 1
			source1 = source[:,::-1]
			source = np.append(source1, source, 1)
			source = np.append(source, source1, 1)
			# axis 2
			source1 = source[:,:,::-1]
			source = np.append(source1, source, 2)
			source = np.append(source, source1, 2)
			# result
			source = source[n0[0]-n1[0]:n1[0]-n0[0],n0[1]-n1[1]:n1[1]-n0[1],n0[2]-n1[2]:n1[2]-n0[2]]
			source1 = 0 #@
			# For beam, always add 0
			beam1 = np.zeros(n0+2*n1, beam.dtype)
			beam1[n0[0]-n1[0]:n1[0]-n0[0],n0[1]-n1[1]:n1[1]-n0[1],n0[2]-n1[2]:n1[2]-n0[2]] = beam
			beam = beam1
			image = np.fft.fftshift(ifft(fft(source)*fft(beam))).real
			image = image[n1[0]:-n1[0],n1[1]:-n1[1],n1[2]:-n1[2]]

	if   (dim == 1) : image = image[:shape[0]]
	elif (dim == 2) : image = image[:shape[0],:shape[1]]
	elif (dim == 3) : image = image[:shape[0],:shape[1],:shape[2]]
	return image

