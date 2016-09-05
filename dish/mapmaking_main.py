#! /usr/bin/env python
import jizhipy as jp
import numpy as np
from jizhipy.Plot import *

from AntArray import *
from Masking import *
#from CaliPhase import *
#from CaliTemperature import *
#from MapMaking import *
##################################################


print 'Start:', jp.Time(1)
outdir = jp.Outdir(0, 'file')


#hdf5dir = '/disk/disk6/SunCygnusA_20160603190723_20160606080722/'
#pixstart, pixlength, pixperiod = 2324, 20, 2564-2324
#sourcename = 'CygA'
#brightsource = jp.BrightSource()


#hdf5dir = '~/study/lal/paon4/paon4_data/CygA/CygA665S1dec15_all-1_ac_0-4095-8'
hdf5dir = '/project/huangqizhi/paon4_data/CygA/CygA665S1dec15_all-1_ac_0-4095-8'
sourcename = 'CygA'
brightsource = jp.BrightSource()


##################################################


antarray = AntArray(hdf5dir)
antarray.WhichHdf5('transitsource')

antarray.SelectVisType('cross1')
#antarray.SelectChannel([1,3,5,7,11,13,15,17,19,25,29,31])
#antarray.SelectChannel([3,5])
#antarray.MaskChannel([9,10,21,23,24,27,28])


##################################################


masking = Masking(antarray)
try : masking.MaskNoiseSource(pixstart, pixlength, pixperiod)
except : pass


freq = 1400
nfreq = abs(antarray.Ant.freq-freq)
nfreq = np.where(nfreq==nfreq.min())[0][0]


masking.See(nfreq, None, 60, 1, 10000, None, 7, 1, True)
masking.See(nfreq, (4000,10000), 60, 1, 10000, None, 7, 1, True)
jp.Raise()


#vis1400 = antarray.vis[:,nfreq,antarray.visorder].flatten()
#np.save('vis1400_data', vis1400)

#nsigma, nloop, threshold, multipool = 5, 20, 0.001, False
#nsigma, nloop, threshold, multipool = 5, None, 0.001, False
#nsigma, nloop, threshold, multipool = 4, None, 0.001, True
#nsigma, nloop, threshold, multipool = 4, 9, 0.001, True
#nsigma, nloop, threshold, multipool = 4, 8, 0.001, True
#nsigma, nloop, threshold, multipool = 4, 7, 0.001, True
#nsigma, nloop, threshold, multipool = 4, 6, 0.001, True
#nsigma, nloop, threshold, multipool = 4, 5, 0.001, True
#nsigma, nloop, threshold, multipool = 4, 4, 0.001, True
#nsigma, nloop, threshold, multipool = 4, 3, 0.001, True
#nsigma, nloop, threshold, multipool = 4, 2, 0.001, True
nsigma, nloop, threshold, multipool = 4, 1, 0.001, False
print 'nsigma='+str(nsigma)+'   nloop='+str(nloop)+'   threshold='+str(threshold)+'   multipool='+str(multipool)
masking.MaskLoop(nsigma=nsigma, nloop=nloop, threshold=threshold, multipool=multipool)

maskmax = masking.mask.sum(0)[:,0]
nfreq_maskmax = np.where(maskmax==maskmax.max())[0][0]
print 'nfreq_maskmax =', nfreq_maskmax
#mask = masking.mask[:,nfreq].flatten()

vis1400 = antarray.vis[:,nfreq_maskmax,antarray.visorder].flatten()
mask = masking.mask[:,nfreq_maskmax].flatten()

#np.save('vis1400_mask_'+str(nsigma)+'-'+str(nloop)+'-'+str(threshold), mask)

vismask = np.ma.MaskedArray(vis1400, mask)

x = np.arange(vis1400.size)
plt.plot(x, vis1400.real, 'b-')
plt.plot(x, vismask.real, 'r-')
plt.savefig('vis1400-max_'+str(nsigma)+'-'+str(nloop)+'-'+str(threshold)+'.png')
exit()


##################################################

caliphase = CaliPhase(antarray, masking)
RA, Dec = brightsource.RADec(sourcename)
print sourcename, '  RA =', RA, '  Dec =', Dec
caliphase.RADec(RA, Dec)
caliphase.Fringe(6)

#caliphase.Smooth(20, 3)
caliphase.FitBeam()
caliphase.FitVis()
caliphase.Plot(1)
Raise()

##################################################

# nsource = 9
#calitemp = CaliTemperature(antarray, 750, (1,25))
calitemp = CaliTemperature(antarray, 750, (8,12))
calitemp.Vis()

pixstart, pixlength, pixperiod = masking.noisesource.pixstart, masking.noisesource.pixlength, masking.noisesource.pixperiod
calitemp.MaskNoiseSource(pixstart, pixlength, pixperiod)

maskmanual = np.zeros(calitemp.vis.shape, bool)
maskmanual[1500:] = True
calitemp.MaskManual(maskmanual)
maskmanual = 0 #@

calitemp.Gaint()
calitemp.vis.mask -= calitemp.masking.maskmanual
calitemp.vis.mask += calitemp.masking.masknoisesource

flux = brightsource.FluxDensity(sourcename, calitemp.freq)

calitemp.Smooth(200)
calitemp.CaliTemp(flux)
Raise()

calitemp.Smooth()
#calitemp.Plot()
calitemp.Smooth(200)

##################################################

mapmaking = MapMaking(caliphase, calitemp)
mapmaking.CaliPhase()
mapmaking.MapMaking(201)


print 'End:', jp.Time(1)
