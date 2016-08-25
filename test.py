#! /usr/bin/env python
import jizhipy as jp

for i in xrange(3) : 

	def check() :

		def moref() : 
			outdir = jp.Outdir(9, 'file')
			return outdir

		outdir = moref()
		return outdir

	outdir = check()

print outdir
