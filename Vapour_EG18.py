#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 15 14:53:14 2019

@author: swa048
"""

import numpy as np
import pandas as pd
from scipy import signal
from scipy import stats as scistats
import datetime as dt
import matplotlib.pyplot as plt

from loadandconvert import pickperiod,loadperiod,loadfluxperiod,removeoutliers,dfperiodfun

pathplot='/Users/swa048/Documents/OneDrive/Dokumente/EGRIP/EGRIP2018/EGRIP18_plots/'
# import the Vapour data

iso7m=pd.read_csv('/Volumes/gfi_snowiso/Vapour/2018/processed_data_combined/EG18_level1.txt',index_col=0,parse_dates=True,na_values=['NAN'])  #7m?
iso30=pd.read_csv('/Volumes/gfi_snowiso/Vapour/2018/processed_data_combined/EG18_level2.txt',index_col=0,parse_dates=True,na_values=['NAN'])  
iso80=pd.read_csv('/Volumes/gfi_snowiso/Vapour/2018/processed_data_combined/EG18_level3.txt',index_col=0,parse_dates=True,na_values=['NAN'])
iso180=pd.read_csv('/Volumes/gfi_snowiso/Vapour/2018/processed_data_combined/EG18_level4.txt',index_col=0,parse_dates=True,na_values=['NAN'])

#import meterological data
''' pick 10 min meterological data for Clear Diurnal Cycle (CDC) period 
     DOY 140-145   Date: 20.05.18-25.05.18 '''
     
CDC_temp10min=pickperiod('/Users/swa048/forServer/Meteo/Ambient/2018/EGRIP/data/EGRIP_TC/','TC_EGRIP18_calibrated.txt','2018-05-20 00:00','2018-05-26 00:00') #calibrated data
CDC_wind10min=loadperiod('/Users/swa048/forServer/Meteo/Ambient/2018/EGRIP/data/','CR3000_Table_ae.dat','2018-05-20 00:00','2018-05-26 00:00')
CDC_wind10min.rename(columns={'WS_ms_Avg(1)': 'wind1_mean','WS_ms_Avg(2)': 'wind2_mean','WS_ms_Avg(3)': 'wind3_mean'}, inplace=True)
CDC_wind10min=CDC_wind10min[['wind1_mean','wind2_mean','wind3_mean']] #wind2 is difficult to interpret



''' pick TK3 processed flux data for Clear Diurnal Cycle (CDC) period 
     DOY 140-145   Date: 20.05.18-25.05.18 '''
CDC_flux=loadfluxperiod('/Users/swa048/forServer/Meteo/EC/2018/EGRIP/fluxes_processed_TK3/','EGRIP2018_period1_result_1901180949.csv','2018-05-20 00:00','2018-05-26 00:00')



''' load the online computed 10 min data'''
fluxON=pd.read_csv('/Users/swa048/forServer/Meteo/EC/2018/EGRIP/onlinefluxes10min/CR3000_ec_scfd.dat',index_col=0,usecols=([0,1,2,3,4,5,24,25,26,27,28]),skiprows=[0,2,3],parse_dates=True,na_values=['NAN',-7999],dtype={'RECORD': np.int64, 'LE': np.float64, 'Hs': np.float64, 'tau': np.float64, 'u_star': np.float64, 'Ux_Ux': np.float64, 'Ux_Uy': np.float64, 'Ux_Uz': np.float64, 'Ux_ln_vh': np.float64, 'Ux_Ts': np.float64, 'Uy_Uy': np.float64, 'Uy_Uz': np.float64, 'Uy_ln_vh': np.float64, 'Uy_Ts': np.float64, 'Uz_Uz': np.float64, 'Uz_ln_vh': np.float64, 'Uz_Ts': np.float64, 'ln_vh_ln_vh': np.float64, 'ln_vh_Ts': np.float64, 'Ts_Ts': np.float64, 'Ux_Avg': np.float64, 'Uy_Avg': np.float64, 'Uz_Avg': np.float64, 'ln_vh_Avg': np.float64, 'Ts_Avg': np.float64, 'horiz_wind_spd': np.float64, 'result_wind_spd': np.float64, 'wind_dir': np.float64, 'wind_dir_std_dev': np.float64, 'wnd_dir_compass': np.float64, 'vh_Avg': np.float64, 'n_Tot': np.float64, 'del_T_f_Tot': np.float64, 'track_f_Tot': np.float64, 'amp_h_f_Tot': np.float64, 'amp_l_f_Tot': np.float64} )
fluxON.rename(columns={'RECORD': 'DOY'}, inplace=True)
t0=np.datetime64(str(fluxON.index[0])[0:4]+'-01-01')
fluxON.loc[:,'DOY'] = pd.Series((fluxON.index-t0)/dt.timedelta(days=1), index=fluxON.index) 

removeoutliers(fluxON,'LE',30,-30)
removeoutliers(iso7m,'d18O',-20,-80)
''' import surface sample'''
pathiso='/Users/swa048/forServer/Snow/Isotopes/2018/data_measuredsamples/ST/'
ST05=pd.read_csv(pathiso+'ST_05cm.txt',delim_whitespace=True,index_col=(0),parse_dates=True)
ST05.DOY=ST05.DOY-1
ST05.index=ST05.index+dt.timedelta(hours=9)
ST1=pd.read_csv(pathiso+'ST_1cm.txt',delim_whitespace=True,index_col=(0),parse_dates=True)
ST1.DOY=ST1.DOY-1
ST1.index=ST1.index+dt.timedelta(hours=9)
ST2=pd.read_csv(pathiso+'ST_2cm.txt',delim_whitespace=True,index_col=(0),parse_dates=True)
ST2.DOY=ST2.DOY-1
ST2.index=ST2.index+dt.timedelta(hours=9)

#add time of precipitation samples
PRECIP=pd.read_csv('/Users/swa048/forServer/Snow/Isotopes/2018/data_measuredsamples/PRECIP/PRECIP_sorted.txt',delim_whitespace=True,index_col=(0),parse_dates=True)
PRECIP.DOY=PRECIP.DOY-1

#
#%% overview 
plt.close('all')

fig_ov,(ax1,ax2)=plt.subplots(2,1,sharex=True,figsize=(16,8))
ax1.plot(iso30.d18O,label='30cm')
ax1.plot(iso80.d18O,label='80cm')
ax1.plot(iso180.d18O,label='180cm')
ax1.plot(iso7m.d18O,label='7m')

ax1.set_ylabel(u'vapour $\delta^{18}$O ‰')

ax1t=ax1.twinx()
ax1t.plot()
ax1t.plot(ST05.d18O,label='0.5cm con',Marker='*',LineStyle='-.',color='violet')
ax1t.plot(ST1.d18O,label='1cm con',Marker='*',LineStyle='-.',color='darkviolet')
ax1t.plot(ST2.d18O,label='2cm con',Marker='*',LineStyle='-.',color='hotpink')
ax1t.plot(PRECIP.d18O,label='precip',Marker='v',Linestyle='',color='darkorange')
ax1t.set_ylabel(u'snow $\delta^{18}$O ‰')
ax1.set_title('EGRIP season 2018')

#ax1.set_ylim(-65,-30)
#ax1t.set_ylim(-37,-27)

ax2.plot(fluxON.LE,color='lightcoral')
ax2.axhline(color='k',linestyle=':')
ax2.set_ylabel('W/m2')
ax2t=ax2.twinx()
ax2t.plot(np.cumsum(fluxON.LE),color='maroon',label='cum LE')
ax2t.set_ylabel('cumulative W/m2')
fig_ov.legend(loc='upper left')
fig_ov.savefig(pathplot+'fullseason.pdf')
#%% clear diurnal cycle  22-07-29-07   203-210
fig_cdc,(ax1,ax2)=plt.subplots(2,1,sharex=True,figsize=(10,8))
ax1.plot(iso30.d18O,label='30cm')
ax1.plot(iso80.d18O,label='80cm')
ax1.plot(iso180.d18O,label='180cm')
ax1.plot(iso7m.d18O,label='7m')
ax1.xaxis.grid(True)
ax1t=ax1.twinx()
ax1t.plot()
ax1t.plot(ST05.d18O,label='0.5cm con',Marker='*',LineStyle='-.',color='violet')
ax1t.plot(ST1.d18O,label='1cm con',Marker='*',LineStyle='-.',color='darkviolet')
ax1t.plot(ST2.d18O,label='2cm con',Marker='*',LineStyle='-.',color='hotpink')
ax1t.plot(PRECIP.d18O,label='precip',Marker='v',Linestyle='',color='darkorange')
ax1.set_ylim(-55,-40)
ax1t.set_ylim(-35,-30)
ax1.set_ylabel(u'vapour $\delta^{18}$O ‰')
ax1t.set_ylabel(u'snow $\delta^{18}$O ‰')


ax2.plot(fluxON.LE,color='lightcoral')
ax2.set_ylabel('W/m2')
ax2.axhline(color='k',linestyle=':')
ax2t=ax2.twinx()
ax2t.plot(np.cumsum(fluxON.LE['2018-07-22':'2018-07-29']),color='maroon',label='cum LE')
ax2t.set_ylabel('cumulative W/m2')
ax1.set_xlim('2018-07-22','2018-07-29')
ax2.set_xlim('2018-07-22','2018-07-29')
fig_cdc.legend(loc='lower center')
fig_cdc.tight_layout()
ax1.set_title('DOY 203-210, no precip')
fig_cdc.savefig(pathplot+'doy203210.pdf')
#%%
hour7=iso7m.d18O.resample('1H').mean()
hour30=iso30.d18O.resample('1H').mean()
hour80=iso80.d18O.resample('1H').mean()
hour180=iso180.d18O.resample('1H').mean()


fig_cdc2,(ax1,ax2)=plt.subplots(2,1,sharex=True,figsize=(10,8))
ax1.plot(hour7-hour30,label='7m-30cm')
ax1.plot(hour7-hour80,label='7m-80cm')
ax1.plot(hour7-hour180,label='7m-180cm')
ax1.axhline(color='k',linestyle=':')
ax1.set_ylabel(u'vapour $\delta^{18}$O ‰')
#ax1.plot(iso7m.d18O,label='7m')
ax1t=ax1.twinx()
ax1t.plot()
ax1t.plot(ST05.d18O,label='0.5cm con',Marker='*',LineStyle='-.',color='violet')
ax1t.plot(ST1.d18O,label='1cm con',Marker='*',LineStyle='-.',color='darkviolet')
ax1t.plot(ST2.d18O,label='2cm con',Marker='*',LineStyle='-.',color='hotpink')
ax1t.plot(PRECIP.d18O,label='precip',Marker='v',Linestyle='',color='darkorange')
ax1t.set_ylabel(u'snow $\delta^{18}$O ‰')
ax1.set_ylim(-3,3)
ax1t.set_ylim(-32,-27)
ax1.text(0.15,0.1,'neg gradient with height', horizontalalignment='center',verticalalignment='center', transform=ax1.transAxes)
ax1.text(0.15,0.9,'pos gradient with height', horizontalalignment='center',verticalalignment='center', transform=ax1.transAxes)
ax1.xaxis.grid(True)

ax2.plot(fluxON.LE,color='lightcoral')
ax2.set_ylabel('W/m2')
ax2t=ax2.twinx()
ax2t.plot(np.cumsum(fluxON.LE['2018-07-22':'2018-07-29']),color='maroon',label='cum LE')
ax2.axhline(color='k',linestyle=':')
ax2t.set_ylabel('cumulative W/m2')
ax1.set_xlim('2018-07-22','2018-07-29')
ax2.set_xlim('2018-07-22','2018-07-29')
fig_cdc2.legend(loc='lower center')
fig_cdc2.tight_layout()
ax1.set_title('DOY 203-210, no precip, vapour gradient')
ax2.xaxis.grid(True)
fig_cdc2.savefig(pathplot+'doy203210_grad.pdf')

#%% overview d-excess


fig_ovex,(ax1,ax2)=plt.subplots(2,1,sharex=True,figsize=(16,8))
ax1.plot(iso30.dD-8*iso30.d18O,label='30cm')
ax1.plot(iso80.dD-8*iso80.d18O,label='80cm')
ax1.plot(iso180.dD-8*iso180.d18O,label='180cm')
ax1.plot(iso7m.dD-8*iso7m.d18O,label='7m')

ax1.set_ylabel(u'vapour d-excess ‰')

ax1t=ax1.twinx()
ax1t.plot()
ax1t.plot(ST05.dexcess,label='0.5cm con',Marker='*',LineStyle='-.',color='violet')
ax1t.plot(ST1.dexcess,label='1cm con',Marker='*',LineStyle='-.',color='darkviolet')
ax1t.plot(ST2.dexcess,label='2cm con',Marker='*',LineStyle='-.',color='hotpink')
ax1t.plot(PRECIP.dexcess,label='precip',Marker='v',Linestyle='',color='darkorange')
ax1t.set_ylabel(u'snow d-excess ‰')
ax1.set_title('d-excess EGRIP season 2018')

#ax1.set_ylim(-65,-30)
#ax1t.set_ylim(-37,-27)

ax2.plot(fluxON.LE,color='lightcoral')
ax2.axhline(color='k',linestyle=':')
ax2.set_ylabel('W/m2')
ax2t=ax2.twinx()
ax2t.plot(np.cumsum(fluxON.LE),color='maroon',label='cum LE')
ax2t.set_ylabel('cumulative W/m2')
fig_ovex.legend(loc='upper left')
fig_ovex.savefig(pathplot+'fullseason_dex.pdf')
#%% clear diurnal cycle  22-07-29-07   203-210
fig_cdcex,(ax1,ax2)=plt.subplots(2,1,sharex=True,figsize=(10,8))
ax1.plot(iso30.dD-8*iso30.d18O,label='30cm')
ax1.plot(iso80.dD-8*iso80.d18O,label='80cm')
ax1.plot(iso180.dD-8*iso180.d18O,label='180cm')
ax1.plot(iso7m.dD-8*iso7m.d18O,label='7m')
ax1.xaxis.grid(True)
ax1t=ax1.twinx()
ax1t.plot()
ax1t.plot(ST05.dexcess,label='0.5cm con',Marker='*',LineStyle='-.',color='violet')
ax1t.plot(ST1.dexcess,label='1cm con',Marker='*',LineStyle='-.',color='darkviolet')
ax1t.plot(ST2.dexcess,label='2cm con',Marker='*',LineStyle='-.',color='hotpink')
ax1t.plot(PRECIP.dexcess,label='precip',Marker='v',Linestyle='',color='darkorange')
ax1.set_ylim(20,45)
ax1t.set_ylim(5,15)
ax1.set_ylabel(u'vapour d-excess ‰')
ax1t.set_ylabel(u'snow d-excess ‰')


ax2.plot(fluxON.LE,color='lightcoral')
ax2.set_ylabel('W/m2')
ax2.axhline(color='k',linestyle=':')
ax2t=ax2.twinx()
ax2t.plot(np.cumsum(fluxON.LE['2018-07-22':'2018-07-29']),color='maroon',label='cum LE')
ax2t.set_ylabel('cumulative W/m2')
ax1.set_xlim('2018-07-22','2018-07-29')
ax2.set_xlim('2018-07-22','2018-07-29')
fig_cdc.legend(loc='lower center')
fig_cdc.tight_layout()
ax1.set_title('d-excess DOY 203-210, no precip')
fig_cdcex.savefig(pathplot+'doy203210_dex.pdf')
#%%
hour7d=iso7m.dD.resample('1H').mean()-8*iso7m.d18O.resample('1H').mean()
hour30d=iso30.dD.resample('1H').mean()-8*iso30.d18O.resample('1H').mean()
hour80d=iso80.dD.resample('1H').mean()-8*iso80.d18O.resample('1H').mean()
hour180d=iso180.dD.resample('1H').mean()-8*iso180.d18O.resample('1H').mean()


fig_cdc2ex,(ax1,ax2)=plt.subplots(2,1,sharex=True,figsize=(10,8))
ax1.plot(hour7d-hour30d,label='7m-30cm')
ax1.plot(hour7d-hour80d,label='7m-80cm')
ax1.plot(hour7d-hour180d,label='7m-180cm')
ax1.axhline(color='k',linestyle=':')
ax1.set_ylabel(u'vapour d-excess ‰')
#ax1.plot(iso7m.d18O,label='7m')
ax1t=ax1.twinx()
ax1t.plot()
ax1t.plot(ST05.dexcess,label='0.5cm con',Marker='*',LineStyle='-.',color='violet')
ax1t.plot(ST1.dexcess,label='1cm con',Marker='*',LineStyle='-.',color='darkviolet')
ax1t.plot(ST2.dexcess,label='2cm con',Marker='*',LineStyle='-.',color='hotpink')
ax1t.plot(PRECIP.dexcess,label='precip',Marker='v',Linestyle='',color='darkorange')
ax1t.set_ylabel(u'snow d-excess ‰')
ax1.set_ylim(-5,10)
ax1t.set_ylim(0,15)
ax1.text(0.15,0.1,'neg gradient with height', horizontalalignment='center',verticalalignment='center', transform=ax1.transAxes)
ax1.text(0.15,0.9,'pos gradient with height', horizontalalignment='center',verticalalignment='center', transform=ax1.transAxes)
ax1.xaxis.grid(True)

ax2.plot(fluxON.LE,color='lightcoral')
ax2.set_ylabel('W/m2')
ax2t=ax2.twinx()
ax2t.plot(np.cumsum(fluxON.LE['2018-07-22':'2018-07-29']),color='maroon',label='cum LE')
ax2.axhline(color='k',linestyle=':')
ax2t.set_ylabel('cumulative W/m2')
ax1.set_xlim('2018-07-22','2018-07-29')
ax2.set_xlim('2018-07-22','2018-07-29')
fig_cdc2.legend(loc='lower center')
fig_cdc2.tight_layout()
ax1.set_title('d-excess DOY 203-210, no precip, vapour gradient')
ax2.xaxis.grid(True)
fig_cdc2ex.savefig(pathplot+'doy203210_graddex.pdf')