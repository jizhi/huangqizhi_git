#! /usr/bin/env python
import jizhipy as jp
import numpy as np
from jizhipy.Plot import *

from AntArray import *
from Masking import *
from CaliGain import *
from CaliPhase import *
#from CaliTemperature import *
#from MapMaking import *
##################################################


print 'Start:', jp.Time(1)+'\n'
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
#antarray.SelectChannel([3, 7])
#antarray.SelectChannel([1, 3, 7])

#antarray.SelectVisType('auto1')
#antarray.SelectChannel([3])


##################################################


masking = Masking(antarray)


freq = 1400
nfreq = abs(antarray.Ant.freq-freq)
nfreq = np.where(nfreq==nfreq.min())[0][0]  # 307
ntimemax = 69104


axis, per, times = 0, 60, 1
#axis, per, times = 1, 7, 1
nsigma, nloop, threshold = 4, 2, 0.001

#masking.MaskLoop(axis=axis, per=per, times=times, nsigma=nsigma, nloop=nloop, threshold=threshold)


#vis = antarray.vis[:,:,antarray.visorder][:,:,None]
##vis[masking.masknoisesource] += 3e4
#x = np.arange(len(vis))
#plt.plot(x, vis[:,nfreq,0], 'b-')
#vis[masking.mask] = masking.maskvalue
#plt.plot(x, vis[:,nfreq,0], 'r-')
#plt.show()
#jp.Raise()


#vis = antarray.vis[:,:,antarray.visorder][:,:,None]
##vis[masking.masknoisesource] += 3e4
#x = np.arange(vis.shape[1])
#plt.plot(x, vis[40000,:,0], 'b-')
#vis[masking.mask] = masking.maskvalue
#plt.plot(x, vis[40000,:,0], 'r-')
#plt.show()
#jp.Raise()


##################################################


caligain = CaliGain(antarray, masking)

caligain.Window(60*4)

#caligain.Gainnu()  # OK, very good
#jp.Raise()

###caligain.See(3, 100, nfreq)

masktime = (5000,11000)
gainttimes = 100

#caligain.Gaint(masktime=masktime, gainttimes=gainttimes, legendsize=4, legendloc=8)


##################################################


nsigma = 4
fwhm2sigmafactor = 2.287

#caliphase = CaliPhase(antarray, masking, caligain, nsigma=nsigma, fwhm2sigmafactor=fwhm2sigmafactor, Nprocess=200)
caliphase = CaliPhase(antarray, None, None, nsigma=nsigma, fwhm2sigmafactor=fwhm2sigmafactor, Nprocess=200, plotfreq=1400)

RA, Dec = brightsource.RADec(sourcename)
print sourcename, '  RA =', RA, '  Dec =', Dec
caliphase.RADec(RA, Dec)


#for nsigma in [2, 3, 4, 5] : 
for nsigma in [4] : 
	caliphase.nsigma = nsigma
	caliphase.Fringe()
	caliphase.Smooth(200)

	
	#vis = caliphase.vis
	#np.save('vis.data', vis.data)
	#np.save('vis.mask', vis.mask)
	#np.save('timem', caliphase.timem)
	#jp.Raise()
	
	#vis = np.load('vis.data.npy')
	#mask = np.load('vis.mask.npy')
	#timem = np.load('timem.npy')
	#vis = np.ma.MaskedArray(vis, mask)
	#caliphase.vis = vis
	#caliphase.timem = timem
	
	
	caliphase.FitBeam()
	#caliphase.Ampb = np.load('mapmaking_main_output/CaliPhase_output/Ampb.npy')
	#caliphase.Timeb = np.load('mapmaking_main_output/CaliPhase_output/Timeb.npy')
	#caliphase.Sigmab = np.load('mapmaking_main_output/CaliPhase_output/Sigmab.npy')
	#caliphase.Phasens = np.load('mapmaking_main_output/CaliPhase_output/Phasens.npy')
	
	
	#caliphase.Nprocess = 2
	caliphase.FitPhase()
	#jp.Raise()
	#caliphase.Lewp = np.load('mapmaking_main_output/CaliPhase_output/Lewp.npy')
	#caliphase.Phaseaddp = np.load('mapmaking_main_output/CaliPhase_output/Phaseaddp_4.npy')
	
	
	#caliphase.FitVis()
	
	caliphase.Plot(dyDeff=0.5, dyLew=0.5, Nprocess=1)


##################################################


# nsource = 9
#calitemp = CaliTemperature(antarray, 750, (1,25))
#calitemp = CaliTemperature(antarray, 750, (8,12))
calitemp = CaliTemperature(antarray, antarray.Ant.freq.mean(), None)
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

calitemp.Smooth(201)
calitemp.CaliTemp(flux)
Raise()

calitemp.Smooth()
#calitemp.Plot()
calitemp.Smooth(201)

##################################################

mapmaking = MapMaking(caliphase, calitemp)
mapmaking.CaliPhase()
mapmaking.MapMaking(201)


print 'End:', jp.Time(1)
