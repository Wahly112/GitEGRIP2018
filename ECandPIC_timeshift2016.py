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

from loadandconvert import loadEC, loadPICfileCH4,loadperiod

#%% compare Tsonic and PROMICE temperature to see if timeseries there is in accordance with each other

#PROMICE=loadperiod('/Users/swa048/forServer/Meteo/PROMICE:GEUS/2016/data/','TOA5_15213.TableMem.dat','2016-06-15 00:00','2016-07-19 12:00')
PROMICE=pd.read_csv('/Users/swa048/forServer/Meteo/PROMICE:GEUS/2016/data/TOA5_15213.TableMem.dat',index_col=0,parse_dates=True,skiprows=[1,2],na_values=['NAN'])

start=PROMICE.index>='2016-06-15 00:00'
end=PROMICE.index<='2016-07-19 12:00'
period=[all(tup) for tup in zip(start,end)]
PROMICE=PROMICE[period]

PROMICE.rename(columns={'Temperature_Avg': 'Tpromice'}, inplace=True)

EC_EGRIP=loadperiod('/Users/swa048/forServer/Meteo/EC/2016/data/','CR3000_ec_scfd.dat','2016-06-15 00:00','2018-08-03 23:59')
EC_EGRIP.rename(columns={'Ts_Avg': 'Tsonic'}, inplace=True)
index_list= EC_EGRIP[EC_EGRIP.ln_vh_Avg <= 6.25].index.tolist()
EC_EGRIP=EC_EGRIP.drop(index_list)
EC_EGRIP['vh_Avg']=np.exp(EC_EGRIP.ln_vh_Avg)
#%%
figtemp,ax1=plt.subplots()
l1=ax1.plot(EC_EGRIP.index,EC_EGRIP.Tsonic,color='red',label='T SONIC 10 min average')
l2=ax1.plot(PROMICE.index,PROMICE.Tpromice,color='blue',label='T PROMICE 10 min average')
ax2=ax1.twinx()
l3=ax2.plot(EC_EGRIP.index,EC_EGRIP.vh_Avg,color='orange',label='KH20 10 min average ')
legs=l1+l2+l3
labs=[l.get_label() for l in legs]
ax1.legend(legs,labs,loc=0)
ax1.set_title('EGRIP 2016')
#ax1.set_xlim('2016-06-18','2016-06-21')

#figtemp.savefig('/Users/swa048/Documents/OneDrive/Dokumente/EC/Isoflux/Isoflux_plots/'+'COMPTsonicTpromice'+'.pdf')






#%% direct comparison
dfPIC1806,fa1806=loadPICfileCH4('/Users/swa048/forServer/Vapour/2016/data/06/18/HIDS2143-20160618-000007Z-DataLog_User.dat','00:00','23:59')
dfPIC1906,fa1906=loadPICfileCH4('/Users/swa048/forServer/Vapour/2016/data/06/19/HIDS2143-20160619-000006Z-DataLog_User.dat','00:00','23:59')
dfPIC2006,fa2006=loadPICfileCH4('/Users/swa048/forServer/Vapour/2016/data/06/20/HIDS2143-20160620-000007Z-DataLog_User.dat','00:00','23:59')

dfPIC=pd.concat((dfPIC1806,dfPIC1906,dfPIC2006),axis=0)

##cut calibration data
#dfPIC=dfPICall[dfPICall.Valve==8]
#index_list= dfPIC[(dfPIC.index >= "2018-06-18 12:30:07") & (dfPIC.index <= "2018-06-18 12:30:30")].index.tolist()
#dfPIC=dfPIC.drop(index_list)
##idxdrop2=pd.date_range(start='2018-06-18 12:30:07',end='2018-06-18 12:30:30',freq='1S')   #manually remove values that stem from calibration remainer
#index_list= dfPIC[(dfPIC.index >= "2018-06-20 16:23:46") & (dfPIC.index <= "2018-06-20 16:24:07")].index.tolist()
#dfPIC=dfPIC.drop(index_list)
##idxdrop=pd.date_range(start='2018-06-20 16:23:46',end='2018-06-20 16:24:07',freq='1S')   #manually remove values that stem from calibration remainer

index_list= dfPIC[(dfPIC.q >= 6000) | (dfPIC.q <= 1600)].index.tolist()
dfPIC=dfPIC.drop(index_list)


#dfPIC=dfPIC.drop(idxdrop)
#dfPIC=dfPIC.drop(idxdrop2)

#shift PIC(UTC) to Local time (-2 hours)
#dfPIC=dfPIC.shift(-2,freq='H')

#
#dfEC1806=loadEC('/Users/swa048/forServer/Meteo/EC/2018/EGRIP/raw20Hzdata/EGRIP_EC_2018169.dat')  #!!! watch out... DOY in datafile is Julian Days  ... should be 168 in filename
#dfEC1906=loadEC('/Users/swa048/forServer/Meteo/EC/2018/EGRIP/raw20Hzdata/EGRIP_EC_2018170.dat')
#dfEC2006=loadEC('/Users/swa048/forServer/Meteo/EC/2018/EGRIP/raw20Hzdata/EGRIP_EC_2018171.dat')
#
#dfEC=pd.concat((dfEC1806,dfEC1906,dfEC2006),axis=0)
#
##resampling the ECq data into 1s means
#ECq1smean=dfEC.q.resample('1s').mean().rename('kh20q')
#ECtemp1smean=dfEC.Ts.resample('1s').mean()

#resampling the ECq data into 1s nearest neighbour
#ECq1s=dfEC.q.resample('1s').first().rename('kh20q')
##%%
#fig1,(ax1,ax2)=plt.subplots(2,1)
#l1=ax1.plot(dfEC.index,dfEC.q,'o',Markersize=1,label='KH20 raw')
##ax2=ax1.twinx()
#l3=ax2.plot(dfPIC.index,dfPIC.q,'o',Markersize=2,label='PIC',color='red')
##l4=ax2.plot(dfPIC_raw.q,label='PIC raw',color='green')
##legs=l1+l2+l3
##labs=[l.get_label() for l in legs]
##ax1.legend(legs,labs,loc=0)
#ax1.legend()
#ax2.legend()
##ax1.xaxis.set_major_formatter(timefmt)
#ax1.set_xlabel('local time')
#ax1.set_title('EGIRP, copper tube '+str(ECq1smean.index[0])[:10])
#limits=ax1.get_xlim()
#ax2.set_xlim(limits)
##fig1.savefig('/Users/swa048/Documents/OneDrive/Dokumente/EC/Isoflux/Isoflux_plots/'+'hum_comparison_EGRIP_1'+'.pdf')
#
#fig2,(ax1,ax3)=plt.subplots(2,1)
#l1=ax1.plot(ECq1smean,'o',Markersize=1,label='KH20 q 1s mean',color='orange')
#ax2=ax1.twinx()
#l3=ax2.plot(dfPIC.index,dfPIC.q,'o',Markersize=1,label='PIC hum not interpolated ',color='blue')
#l4=ax3.plot(ECtemp1smean,color='red',label='Tsonic')
##l4=ax2.plot(dfPIC_raw.q,label='PIC raw',color='green')
#legs=l1+l4+l3
#labs=[l.get_label() for l in legs]
#ax1.legend(legs,labs,loc=0)
##ax1.legend()
#
##ax1.xaxis.set_major_formatter(timefmt)
#ax1.set_xlabel('local time')
#ax1.set_title('EGIRP, copper tube '+str(ECq1smean.index[0])[:10])
##limits=ax1.get_xlim()
##ax2.set_xlim(limits)
#
#fig3,ax1=plt.subplots()
##l1=ax1.plot(ECq1smean,'o',Markersize=1,label='KH20 q 1s mean',color='orange')
#ax2=ax1.twinx()
#l3=ax2.plot(dfPIC.index,dfPIC.q,'o',Markersize=1,label='PIC hum not interpolated ',color='blue')
#l4=ax1.plot(ECtemp1smean,color='red',label='Tsonic')
##l4=ax2.plot(dfPIC_raw.q,label='PIC raw',color='green')
#legs=l4+l3
#labs=[l.get_label() for l in legs]
#ax1.legend(legs,labs,loc=0)
##ax1.legend()
#
##ax1.xaxis.set_major_formatter(timefmt)
#ax1.set_xlabel('local time')
#ax1.set_title('EGIRP, copper tube '+str(ECq1smean.index[0])[:10])
##limits=ax1.get_xlim()
##ax2.set_xlim(limits)

fig5,ax1=plt.subplots()
l4=ax1.plot(EC_EGRIP.index[(EC_EGRIP.index >= "2016-06-18 00:00") & (EC_EGRIP.index <= "2016-06-21 00:00")],EC_EGRIP.vh_Avg[(EC_EGRIP.index >= "2016-06-18 00:00") & (EC_EGRIP.index <= "2016-06-21 00:00")],label='KH20 10 min average ',color='red')
ax2=ax1.twinx()
l3=ax2.plot(dfPIC.index,dfPIC.q,'o',Markersize=1,label='PIC hum not interpolated ',color='blue')
legs=l4+l3
labs=[l.get_label() for l in legs]
ax1.legend(legs,labs,loc=0)
ax1.set_title('EGRIP 2016')

#ax1.xaxis.set_major_formatter(timefmt)
ax1.set_xlabel('UTC')
