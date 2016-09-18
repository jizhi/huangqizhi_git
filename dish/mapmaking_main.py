#! /usr/bin/env python
import jizhipy as jp
import numpy as np
from jizhipy.Plot import *

from AntArray import *
from Masking import *
from CaliGain import *
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
antarray.SelectChannel([3,5])

#antarray.SelectVisType('auto1')
#antarray.SelectChannel([3])


##################################################


masking = Masking(antarray)
try : masking.MaskNoiseSource(pixstart, pixlength, pixperiod)
except : pass


freq = 1400
nfreq = abs(antarray.Ant.freq-freq)
nfreq = np.where(nfreq==nfreq.min())[0][0]
ntimemax = 69104


timeper, timetimes, freqper, freqtimes = 60, 1, 0, 1
nsigma, nloop, threshold, Nprocess = 4, 1, 0.001, None
filtersize, verbose = None, True

masking.MaskLoop(timeper=timeper, timetimes=timetimes, freqper=freqper, freqtimes=freqtimes, nsigma=nsigma, nloop=nloop, threshold=threshold, Nprocess=Nprocess, filtersize=filtersize, verbose=verbose)


##################################################


caligain = CaliGain(antarray, masking)

caligain.Window(4*60)

#caligain.See(3, 100, nfreq)

caligain.Gaint()

caligain.Smooth([3,1], 100)
jp.Raise()









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
