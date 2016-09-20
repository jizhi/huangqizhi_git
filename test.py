#! /usr/bin/env python
import numpy as np
import jizhipy as jp
import time


def ArrayAxis( array, axis1, axis2, act ) : 
	'''
	act:
		'move' or 'exchange'
		If (act == 'move')     : move     a[axis1] to   b[axis2]
		If (act == 'exchange') : exchange a[axis1] with a[axis2]
	'''
	array = npfmt(array)
	shapeo = array.shape
	if (len(shapeo) <= 1) : return array
	# Check axis
	if (axis1 < 0) : axis1 = len(shapeo) + axis1
	if (axis2 < 0) : axis2 = len(shapeo) + axis2
	if (axis1>=len(shapeo) or axis2>=len(shapeo)) : Raise(Exception, 'axis1='+str(axis1)+', axis2='+str(axis2)+' out of array.shape='+str(shapeo)+'=>'+str(len(shapeo))+'D')
	def ArrayAxis2First( array, axis1 ) : 
		if (axis1 == 0) : return array
		array_new = []
		for i in xrange(shapeo[axis1]) : 
			if (axis1 == 1) : 
				array_new = array_new + [array[:,i]]
			elif (axis1 == 2) : 
				array_new = array_new + [array[:,:,i]]
			elif (axis1 == 3) : 
				array_new = array_new + [array[:,:,:,i]]
			elif (axis1 == 4) : 
				array_new = array_new + [array[:,:,:,:,i]]
			elif (axis1 == 5) : 
				array_new = array_new + [array[:,:,:,:,:,i]]
			else : Raise(Exception, 'this function can just handel 6D array. For >= 7D array, you can modify this function by yourself.')
		return npfmt(array_new)
	#@#@#@#@#@#@
	shapeo = list(shapeo)
	s0 = list(np.arange(len(shapeo)))
	act = act.lower()
	if (act == 'move') : 
		if (axis1 < axis2) : 
			s01, s02, s03 = s0[:axis1], s0[axis1+1:axis2+1], s0[axis2+1:]
			s1 = s01 + s02 + [s0[axis1]] + s03
		elif (axis1 > axis2) : 
			s01, s02, s03 = s0[:axis2], s0[axis2:axis1], s0[axis1+1:]
			s1 = s01 + [s0[axis1]] + s02 + s03
		else : return array
	elif (act == 'exchange') : 
		if (axis1 == axis2) : return array
		s1 = s0[:]
		s01, s02 = s0[axis1], s0[axis2]
		s1[axis1], s1[axis2] = s02, s01
	else : Raise(Exception, 'act='+act+' is not "move" nor "exchange"')
	if (len(shapeo) == 1) : return array
	elif (len(shapeo) == 2) : return array.T
	for i in xrange(len(s0)-1, -1, -1) : 
		while (s0[i] != s1[i]) : 
			array = ArrayAxis2First(array, i)
			s0 = [s0[i]] + s0[:i] + s0[i+1:]
			shapeo = [shapeo[i]] + shapeo[:i] + shapeo[i+1:]
	return array




#def ArrayAxis2( array, axis1, axis2, act='move' ) : 
#	array = npfmt(array)
#	shapeo = array.shape
#	if (len(shapeo) <= 1) : return array
#	if (axis1 < 0) : axis1 = len(shapeo) + axis1
#	if (axis2 < 0) : axis2 = len(shapeo) + axis2
#	if (axis1>=len(shapeo) or axis2>=len(shapeo)) : Raise(Exception, 'axis1='+str(axis1)+', axis2='+str(axis2)+' out of array.shape='+str(shapeo)+'=>'+str(len(shapeo))+'D')
#	#--------------------------------------------------
#
#	if (len(shapeo) == 2) : return array.T
#
#	if (len(shapeo) == 3) : 


a = np.arange(2*3*4*5).reshape(2,3,4,5)

