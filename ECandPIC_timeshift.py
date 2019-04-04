#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 18 08:59:31 2019

Check for timeshift in EC and PIC data of EGRIP

@author: swa048
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from loadandconvert import loadEC, loadPIC17,loadperiod

#%% compare Tsonic and PROMICE temperature to see if timeseries there is in accordance with each other

PROMICE=loadperiod('/Users/swa048/forServer/Meteo/PROMICE:GEUS/2018/data_10min_2017+2018/','TOA5_15213.TableMem_3.dat','2018-05-20 00:00','2018-07-29 23:59')
PROMICE.rename(columns={'Temperature_Avg': 'Tpromice'}, inplace=True)

EC_EGRIP=loadperiod('/Users/swa048/forServer/Meteo/EC/2018/EGRIP/onlinefluxes10min/','CR3000_ec_scfd.dat','2018-05-20 00:00','2018-07-29 23:59')
EC_EGRIP.rename(columns={'Ts_Avg': 'Tsonic'}, inplace=True)
index_list= EC_EGRIP[EC_EGRIP.ln_vh_Avg <= 5].index.tolist()
EC_EGRIP=EC_EGRIP.drop(index_list)

figtemp,ax1=plt.subplots()
ax1.plot(EC_EGRIP.index,EC_EGRIP.Tsonic,color='red',label='T SONIC 10 min average')
ax1.plot(PROMICE.index,PROMICE.Tpromice,color='blue',label='T PROMICE 10 min average')
ax2=ax1.twinx()
ax2.plot(EC_EGRIP.index,EC_EGRIP.vh_Avg,color='orange',label='KH20 10 min average ')
ax1.legend()
ax1.set_xlim('2018-06-18','2018-06-21')
#figtemp.savefig('/Users/swa048/Documents/OneDrive/Dokumente/EC/Isoflux/Isoflux_plots/'+'COMPTsonicTpromice'+'.pdf')






#%% direct comparison
dfPIC1806,fa1806=loadPIC17('/Users/swa048/forServer/Vapour/2018/Picarro_2140/data/06/18/HKDS2032-20180618-000006Z-DataLog_User.dat','00:00','23:59')
dfPIC1906,fa1906=loadPIC17('/Users/swa048/forServer/Vapour/2018/Picarro_2140/data/06/19/HKDS2032-20180619-000008Z-DataLog_User.dat','00:00','23:59')
dfPIC2006,fa2006=loadPIC17('/Users/swa048/forServer/Vapour/2018/Picarro_2140/data/06/20/HKDS2032-20180620-000007Z-DataLog_User.dat','00:00','23:59')

dfPICall=pd.concat((dfPIC1806,dfPIC1906,dfPIC2006),axis=0)

#cut calibration data
dfPIC=dfPICall[dfPICall.Valve==8]
index_list= dfPIC[(dfPIC.index >= "2018-06-18 12:30:07") & (dfPIC.index <= "2018-06-18 12:30:30")].index.tolist()
dfPIC=dfPIC.drop(index_list)
#idxdrop2=pd.date_range(start='2018-06-18 12:30:07',end='2018-06-18 12:30:30',freq='1S')   #manually remove values that stem from calibration remainer
index_list= dfPIC[(dfPIC.index >= "2018-06-20 16:23:46") & (dfPIC.index <= "2018-06-20 16:24:07")].index.tolist()
dfPIC=dfPIC.drop(index_list)
#idxdrop=pd.date_range(start='2018-06-20 16:23:46',end='2018-06-20 16:24:07',freq='1S')   #manually remove values that stem from calibration remainer

index_list= dfPIC[(dfPIC.q >= 5000) & (dfPIC.q <= 2000)].index.tolist()
dfPIC=dfPIC.drop(index_list)


#dfPIC=dfPIC.drop(idxdrop)
#dfPIC=dfPIC.drop(idxdrop2)

#shift PIC(UTC) to Local time (-2 hours)
#dfPIC=dfPIC.shift(-2,freq='H')


dfEC1806=loadEC('/Users/swa048/forServer/Meteo/EC/2018/EGRIP/raw20Hzdata/EGRIP_EC_2018169.dat')  #!!! watch out... DOY in datafile is Julian Days  ... should be 168 in filename
dfEC1906=loadEC('/Users/swa048/forServer/Meteo/EC/2018/EGRIP/raw20Hzdata/EGRIP_EC_2018170.dat')
dfEC2006=loadEC('/Users/swa048/forServer/Meteo/EC/2018/EGRIP/raw20Hzdata/EGRIP_EC_2018171.dat')

dfEC=pd.concat((dfEC1806,dfEC1906,dfEC2006),axis=0)

#resampling the ECq data into 1s means
ECq1smean=dfEC.q.resample('1s').mean().rename('kh20q')
ECtemp1smean=dfEC.Ts.resample('1s').mean()

#resampling the ECq data into 1s nearest neighbour
#ECq1s=dfEC.q.resample('1s').first().rename('kh20q')
#%%
fig1,(ax1,ax2)=plt.subplots(2,1)
l1=ax1.plot(dfEC.index,dfEC.q,'o',Markersize=1,label='KH20 raw')
#ax2=ax1.twinx()
l3=ax2.plot(dfPIC.index,dfPIC.q,'o',Markersize=2,label='PIC',color='red')
#l4=ax2.plot(dfPIC_raw.q,label='PIC raw',color='green')
#legs=l1+l2+l3
#labs=[l.get_label() for l in legs]
#ax1.legend(legs,labs,loc=0)
ax1.legend()
ax2.legend()
#ax1.xaxis.set_major_formatter(timefmt)
ax1.set_xlabel('local time')
ax1.set_title('EGIRP, copper tube '+str(ECq1smean.index[0])[:10])
limits=ax1.get_xlim()
ax2.set_xlim(limits)
#fig1.savefig('/Users/swa048/Documents/OneDrive/Dokumente/EC/Isoflux/Isoflux_plots/'+'hum_comparison_EGRIP_1'+'.pdf')

fig2,(ax1,ax3)=plt.subplots(2,1)
l1=ax1.plot(ECq1smean,'o',Markersize=1,label='KH20 q 1s mean',color='orange')
ax2=ax1.twinx()
l3=ax2.plot(dfPIC.index,dfPIC.q,'o',Markersize=1,label='PIC hum not interpolated ',color='blue')
l4=ax3.plot(ECtemp1smean,color='red',label='Tsonic')
#l4=ax2.plot(dfPIC_raw.q,label='PIC raw',color='green')
legs=l1+l4+l3
labs=[l.get_label() for l in legs]
ax1.legend(legs,labs,loc=0)
#ax1.legend()

#ax1.xaxis.set_major_formatter(timefmt)
ax1.set_xlabel('local time')
ax1.set_title('EGIRP, copper tube '+str(ECq1smean.index[0])[:10])
#limits=ax1.get_xlim()
#ax2.set_xlim(limits)

fig3,ax1=plt.subplots()
#l1=ax1.plot(ECq1smean,'o',Markersize=1,label='KH20 q 1s mean',color='orange')
ax2=ax1.twinx()
l3=ax2.plot(dfPIC.index,dfPIC.q,'o',Markersize=1,label='PIC hum not interpolated ',color='blue')
l4=ax1.plot(ECtemp1smean,color='red',label='Tsonic')
#l4=ax2.plot(dfPIC_raw.q,label='PIC raw',color='green')
legs=l4+l3
labs=[l.get_label() for l in legs]
ax1.legend(legs,labs,loc=0)
#ax1.legend()

#ax1.xaxis.set_major_formatter(timefmt)
ax1.set_xlabel('local time')
ax1.set_title('EGIRP, copper tube '+str(ECq1smean.index[0])[:10])
#limits=ax1.get_xlim()
#ax2.set_xlim(limits)

fig5,ax1=plt.subplots()
l1=ax1.plot(ECq1smean,'o',Markersize=1,label='KH20 q 1s mean',color='orange')
ax2=ax1.twinx()
l3=ax2.plot(dfPIC.index,dfPIC.q,'o',Markersize=1,label='PIC hum not interpolated ',color='blue')
l4=ax1.plot(EC_EGRIP.index[(EC_EGRIP.index >= "2018-06-18 00:00") & (EC_EGRIP.index <= "2018-06-21 00:00")],EC_EGRIP.vh_Avg[(EC_EGRIP.index >= "2018-06-18 00:00") & (EC_EGRIP.index <= "2018-06-21 00:00")],label='KH20 10 min average ',color='red')
#l4=ax2.plot(dfPIC_raw.q,label='PIC raw',color='green')
legs=l1+l4+l3
labs=[l.get_label() for l in legs]
ax1.legend(legs,labs,loc=0)
#ax1.legend()

#ax1.xaxis.set_major_formatter(timefmt)
ax1.set_xlabel('local time')
ax1.set_title('EGIRP, copper tube '+str(ECq1smean.index[0])[:10])
#limits=ax1.get_xlim()
#ax2.set_xlim(limits)