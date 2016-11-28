#! /usr/bin/env python
import numpy as np
import jizhipy as jp
import time

obstimestart = '2016/05/28 17:11:28.35'
sec1970start = 1464426688.35

print obstimestart
print jp.Time(obstimestart, 8, 0)





#try : sec19701 = time.mktime(time.strptime(time1[:n], '%Y/%m/%d %p %I:%M:%S')) + dot
#except : sec19701 = time.mktime(time.strptime(time1[:n], '%Y/%m/%d %H:%M:%S')) + dot
#if (time2 is None) : return sec19701
