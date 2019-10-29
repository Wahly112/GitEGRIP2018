#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 15 15:30:03 2019

import TK3 calculated csv files and transform them in one dataframe

@author: swa048
"""
import math
import pandas as pd
import numpy as np
import datetime as dt
from scipy import interpolate
from scipy import signal
from scipy import stats as scistats
import matplotlib.pyplot as plt

from loadandconvert import loadfluxperiod,removeoutliers,loadCSIfluxperiod,dfperiodfun,loadperiod
import scipy.stats

import windrose
import matplotlib.cm as cm

#plt.close('all')

#%% EGRIP season 20.05.18-03.08.18
#read in the TK3 flux files
EGRIP=0
if EGRIP:
    path='/Users/swa048/forServer/Meteo/EC/2018/EGRIP/fluxes_processed_TK3/'
    file1='EGRIP2018_period1_result_1901180949.csv'
    file2='EGRIP2018_period2_result_1901222249.csv'
    file3='EGRIP2018_period3_result_1901181832.csv'
    file12nd='EGRIP2018_period12nd_result_1908071535.csv'
    file22nd='EGRIP2018_period2new_result_1908110115.csv'
    file32nd='EGRIP2018_period3new_result_1908130218.csv'
    
    
    flux1=loadfluxperiod(path,file1,'2018-05-20 00:15','2018-05-25 16:15')
    flux12nd=loadfluxperiod(path,file12nd,'2018-05-20 00:15','2018-05-25 16:15')
    flux2=loadfluxperiod(path,file2,'2018-05-25 17:15','2018-07-04 19:45')
    flux22nd=loadfluxperiod(path,file22nd,'2018-05-25 17:15','2018-07-04 19:45')
    flux3=loadfluxperiod(path,file3,'2018-07-04 20:15','2018-08-03 11:15')
    flux32nd=loadfluxperiod(path,file32nd,'2018-07-04 20:15','2018-08-03 11:15')
    
    fluxTK3=pd.concat((flux1,flux2,flux3),sort=False)#sort=True)
    fluxTK32nd=pd.concat((flux12nd,flux22nd,flux32nd),sort=False)
    #wind direction is already taking the 242° into account
    fluxTK3per=dfperiodfun(fluxTK3,'2018-06-11','2018-06-20')
    
    
    # read in the online computed fluxes
    #the N_offset for wind_dir_compass was here assumed to be 180 when in fact it was 242 (or -118)
    fluxON=pd.read_csv('/Users/swa048/forServer/Meteo/EC/2018/EGRIP/onlinefluxes10min/CR3000_ec_scfd.dat',index_col=0,usecols=([0,1,2,3,4,5,24,25,26,27,28,30]),skiprows=[0,2,3],parse_dates=True,na_values=['NAN',-7999],dtype={'RECORD': np.int64, 'LE': np.float64, 'Hs': np.float64, 'tau': np.float64, 'u_star': np.float64, 'Ux_Ux': np.float64, 'Ux_Uy': np.float64, 'Ux_Uz': np.float64, 'Ux_ln_vh': np.float64, 'Ux_Ts': np.float64, 'Uy_Uy': np.float64, 'Uy_Uz': np.float64, 'Uy_ln_vh': np.float64, 'Uy_Ts': np.float64, 'Uz_Uz': np.float64, 'Uz_ln_vh': np.float64, 'Uz_Ts': np.float64, 'ln_vh_ln_vh': np.float64, 'ln_vh_Ts': np.float64, 'Ts_Ts': np.float64, 'Ux_Avg': np.float64, 'Uy_Avg': np.float64, 'Uz_Avg': np.float64, 'ln_vh_Avg': np.float64, 'Ts_Avg': np.float64, 'horiz_wind_spd': np.float64, 'result_wind_spd': np.float64, 'wind_dir': np.float64, 'wind_dir_std_dev': np.float64, 'wnd_dir_compass': np.float64, 'vh_Avg': np.float64, 'n_Tot': np.float64, 'del_T_f_Tot': np.float64, 'track_f_Tot': np.float64, 'amp_h_f_Tot': np.float64, 'amp_l_f_Tot': np.float64} )
    fluxON.rename(columns={'RECORD': 'DOY',
                           'wnd_dir_compass':'wind_dir_compass',
                           "horiz_wind_spd":'wind_spd'}, inplace=True)
    t0=np.datetime64(str(fluxON.index[0])[0:4]+'-01-01')
    fluxON.loc[:,'DOY'] = pd.Series((fluxON.index-t0)/dt.timedelta(days=1), index=fluxON.index) 
    
    removeoutliers(fluxON,'LE',50,-50)
    removeoutliers(fluxON,'Hs',50,-50)
    fluxONper=dfperiodfun(fluxON,'2018-06-11','2018-06-20')
    
    period='30'
    #read in the computed wind dir and wind speed from the 20Hz data for the time, we coutinuously ran in 2m height
    #period 
    dfWDwq=pd.read_csv('/Users/swa048/Documents/OneDrive/Dokumente/EC/Isoflux/datatables/dfWDwq_'+period+'_EGRIP.txt',index_col=0,parse_dates=True,na_values=['NAN'])    
    dfWSwq=pd.read_csv('/Users/swa048/Documents/OneDrive/Dokumente/EC/Isoflux/datatables/dfWSwq_'+period+'_EGRIP.txt',index_col=0,parse_dates=True,na_values=['NAN'])  
    
    Noffset=242
    
    allWD_compass=(dfWDwq.iloc[0] +Noffset) % 360
    allWD=dfWDwq.iloc[0] % 360
    allWS=dfWSwq.iloc[0]

    #%% pearson correlation
    #latent heat
    #fluxEGRIP18=pd.concat([fluxON[['LE','Hs']],fluxTK3[['LvE','HTs']],fluxTK32nd],axis=1)
    #flux_cleanLE=fluxcomb[['LE','LvE']].dropna()
    #print (scipy.stats.pearsonr(flux_cleanLE.LE,flux_cleanLE.LvE))
    # sensible heat
    #flux_cleanHS=fluxcomb[['Hs','HTs']].dropna()
    #print (scipy.stats.pearsonr(flux_cleanHS.Hs,flux_cleanHS.HTs))

    
    #comparison daily mean values
    EGRIP18dailyold=fluxTK3.resample('1D').mean()  #
    EGRIP18daily=fluxTK32nd.resample('1D').mean()
    fluxONdaily=fluxON.resample('1D').mean()
    
#%%    windrose comparison
    windroseEGRIP=plt.figure(figsize=(24,6))

    #Subplot 1: 
    ax1 = windroseEGRIP.add_subplot(131,projection='windrose')
    ax1.bar(allWD_compass, allWS, normed=True, opening=1, edgecolor='white',cmap=cm.hot)
    ax1.set_legend()
    ax1.set_title('30min own calculation, wind-dir-compass')
    
    #Subplot 2:
    ax2 = windroseEGRIP.add_subplot(132,projection='windrose')
    ax2.bar(fluxONper.wind_dir_compass+62, fluxONper.wind_spd, normed=True, opening=1, edgecolor='white',cmap=cm.hot)
    ax2.set_legend()
    ax2.set_title('Online 10min values, wind-dir-compass')
    
    #Subplot 3:
    ax3 = windroseEGRIP.add_subplot(133,projection='windrose')
    ax3.bar(fluxTK3per.wind_dir, fluxTK3per.u, normed=True, opening=1, edgecolor='white',cmap=cm.hot)
    ax3.set_legend()
    ax3.set_title('TK3 half hour values, wind-dir-compass')
    windroseEGRIP.savefig('/Users/swa048/Documents/OneDrive/Dokumente/EGRIP/EGRIP2018/EGRIP18_plots/wind/'+'windrose_EGRIP18'+'.pdf') 
    
    #%% compare visually
    
    fig1,(ax1,ax2)=plt.subplots(2,1)
    ax1.plot(fluxTK32nd.index,fluxTK32nd.LvE,color='seagreen',label='TK3 k=-0.2')
    ax1.plot(fluxON.index,fluxON.LE,color='fuchsia',label='online')
    ax1.set_ylabel('Wm-2')
    

    ax1.set_title('latent heat')
    ax1.axhline(linestyle=':')
    
#    ax2.plot(fluxON.LE,fluxTK32nd.LvE,Marker='+',LineStyle='',label='comp')
#    ax2.set_xlabel('online flux')
#    ax2.set_ylabel('TK3 k=-0.2')
    
#    ax3.plot(flux_cleanLE.LE-flux_cleanLE.LvE,label='Online-TK3')
#    ax3.set_ylabel('Wm-2')
#    ax3.text(0.2, 0.1, f'mean diff = {np.round(np.mean(flux_cleanLE.LE-flux_cleanLE.LvE),2)}', horizontalalignment='center',verticalalignment='center', transform=ax3.transAxes)
#    ax3.text(0.6, 0.1, f'median diff = {np.round(np.median(flux_cleanLE.LE-flux_cleanLE.LvE),2)}', horizontalalignment='center',verticalalignment='center', transform=ax3.transAxes)
    fig1.legend()
    fig1.savefig('/Users/swa048/Documents/OneDrive/Dokumente/EGRIP/EGRIP2018/EGRIP18_plots/fluxes/'+'fluxLE_comparison_EGRIP'+'.pdf') 
    
    
    #%% compare visually
    
    fig2,(ax1,ax2,ax3)=plt.subplots(3,1)
    ax1.plot(fluxON.Hs,color='darkorchid',label='online')
    ax1.set_ylabel('Wm-2')
    ax1.set_xlim('2018-05-20','2018-08-03')
    
    #ax2=ax1.twinx()
    ax1.plot(fluxTK32nd.HTs,color='palegreen',label='TK3')
    ax1.set_title('sensible heat')
    ax1.axhline(linestyle=':')

    
#    ax2.plot(fluxON.Hs,fluxTK32nd.HTs,Marker='+',LineStyle='',label='comp')
#    ax2.set_xlabel('online flux')
#    ax2.set_ylabel('TK3 new')

    
#    ax3.plot(fluxEGRIP18.Hs-fluxEGRIP18.HTs,label='Online-TK3')
#    ax3.set_ylabel('Wm-2')    
#    ax3.text(0.2, 0.1, f'mean diff = {np.round(np.mean(fluxEGRIP18.Hs-fluxEGRIP18.HTs),2)}', horizontalalignment='center',verticalalignment='center', transform=ax3.transAxes)
#    ax3.text(0.6, 0.1, f'median diff = {np.round(np.median(fluxEGRIP18.Hs-fluxEGRIP18.HTs),2)}', horizontalalignment='center',verticalalignment='center', transform=ax3.transAxes)
#    ax3.set_xlim('2018-05-20','2018-08-03')    
    
    fig2.legend()
    fig2.savefig('/Users/swa048/Documents/OneDrive/Dokumente/EGRIP/EGRIP2018/EGRIP18_plots/fluxes/'+'fluxH_comparison_EGRIP'+'.pdf') 
    #%% compare two different flux calculations in TK3 and the sensitiveiy of the outcome
    figx=plt.figure()
    ax1=plt.subplot(221)
    ax1.plot(EGRIP18dailyold.index,EGRIP18dailyold.LvE,label='KH20 kw=-0.144',color='brown')
    ax1.plot(EGRIP18daily.index,EGRIP18daily.LvE,label='KH20 kw=-0.2',color='gold')
    ax1.plot(fluxON.index,fluxON.LE,label='online',color='red')
    ax1.set_title('timeseries LE')
    ax3=plt.subplot(222)
    ax3.plot(EGRIP18daily.LvE-EGRIP18dailyold.LvE,label='kw0.2-kw0.144',color='brown')
    ax3.plot(EGRIP18daily.LvE-fluxON.LE,label='kw0.2-Online',color='red')
    ax3.set_title('difference to new calculations')
    
    ax2=plt.subplot(223)
    ax2.plot(EGRIP18daily.index,EGRIP18daily.HTs,label='TK3 H',color='gold')
    ax2.plot(EGRIP18dailyold.index,EGRIP18dailyold.HTs,label='online H',color='brown')
    ax2.set_title('timeseries H')
    ax4=plt.subplot(224)
    ax4.plot(EGRIP18daily.HTs-EGRIP18dailyold.HTs,label='kw0.2-kw0.144')
    ax4.set_title('')
    figx.legend()    

    #%%    daily mean comparison
    figdaily,(ax1,ax2,ax3)=plt.subplots(3,1)
    ax1.plot(fluxONdaily.LE,color='darkorchid',label='online',MarkerSize=4)
    ax1.set_ylabel('Wm-2')
    ax1.set_xlim('2018-05-20','2018-08-03')
    
    #ax2=ax1.twinx()
    ax1.plot(EGRIP18daily.LvE,color='palegreen',label='TK3 k=-0.2',MarkerSize=4)
    ax1.set_title('Latent Heat Daily Means EGRIP 2018')
    ax1.axhline(linestyle=':')
    
    ax1.plot(EGRIP18dailyold.index,EGRIP18dailyold.LvE,color='brown',label='TK3 k=-0.144',MarkerSize=4)

    
    ax2.plot(fluxONdaily.LE-EGRIP18dailyold.LvE,label='online-TK3old')
    ax2.plot(fluxONdaily.LE-EGRIP18daily.LvE,label='online - TK3new')
    ax2.set_ylabel(' Wm-2')
    ax2.text(0.2, 0.1, f'mean diffnew = {np.round(np.mean(fluxONdaily.LE-EGRIP18daily.LvE),2)}', horizontalalignment='center',verticalalignment='center', transform=ax3.transAxes)
    ax2.text(0.6, 0.1, f'variance diff new= {np.round(np.var(fluxONdaily.LE-EGRIP18daily.LvE),2)}', horizontalalignment='center',verticalalignment='center', transform=ax3.transAxes)
    ax2.set_xlim('2018-05-20','2018-08-03')
    #ax2.set_ylim(-5,5)
    ax2.grid(True)
    
    ax3.hist(fluxONdaily['2018-05-20':'2018-08-03'].LE-EGRIP18daily.LvE,40,label='diff')
    figdaily.legend()
    figdaily.savefig('/Users/swa048/Documents/OneDrive/Dokumente/EGRIP/EGRIP2018/EGRIP18_plots/fluxes/'+'fluxLE_dailycomparison_EGRIP'+'.pdf') 
    
    
    #%%H sensible heat
    figdaily,(ax1,ax2,ax3)=plt.subplots(3,1)
    ax1.plot(fluxONdaily.Hs,color='darkorchid',label='online',MarkerSize=4)
    ax1.set_ylabel('Wm-2')
    ax1.set_xlim('2018-05-20','2018-08-03')
    
    #ax2=ax1.twinx()
    ax1.plot(EGRIP18daily.HTs,color='palegreen',label='TK3',MarkerSize=4)
    ax1.set_title('sensible Heat Daily Means EGRIP 2018')
    ax1.axhline(linestyle=':')

    
    ax2.plot(fluxONdaily.Hs-EGRIP18daily.HTs,label='online-TK3')
    ax2.set_ylabel(' Wm-2')
    ax2.text(0.2, 0.1, f'mean diff = {np.round(np.mean(fluxONdaily.Hs-EGRIP18daily.HTs),2)}', horizontalalignment='center',verticalalignment='center', transform=ax3.transAxes)
    ax2.text(0.6, 0.1, f'variance diff = {np.round(np.var(fluxONdaily.Hs-EGRIP18daily.HTs),2)}', horizontalalignment='center',verticalalignment='center', transform=ax3.transAxes)
    ax2.set_xlim('2018-05-20','2018-08-03')
    #ax2.set_ylim(-5,5)
    ax2.grid(True)
    
    ax3.hist(fluxONdaily.Hs['2018-05-20':'2018-08-03']-EGRIP18daily.HTs['2018-05-20':'2018-08-03'],40,label='diff')
    figdaily.legend()
    figdaily.savefig('/Users/swa048/Documents/OneDrive/Dokumente/EGRIP/EGRIP2018/EGRIP18_plots/fluxes/'+'fluxHs_dailycomparison_EGRIP'+'.pdf') 
 
#%%  season 13.12.18-22.01.19
KOHNEN=0
if KOHNEN:
    
    pathKO='/Users/swa048/forServer/Meteo/EC/2018/KOHNEN/TK3/'
    file1KO='KOHNEN2018_result_1904182151.csv'
    file2KO='KOHNEN2019_result_1904211414.csv'
    
    NoffsetKo=65 #degrees from North
    
    flux1KO=loadfluxperiod(pathKO,file1KO,'2018-12-13 12:00','2018-12-31 23:59')
    flux2KO=loadfluxperiod(pathKO,file2KO,'2019-01-01 00:00','2019-01-22 23:59')
    
    
    fluxTK3_KO=pd.concat((flux1KO,flux2KO),sort=False)#sort=True)
    
    fluxIRGA=loadCSIfluxperiod('/Users/swa048/forServer/Meteo/EC/2018/KOHNEN/CR6Series_Flux_CSIFormat.dat','2018-12-13 12:00','2019-01-22 23:59')
    
    flux1min=loadperiod('/Users/swa048/forServer/Meteo/EC/2018/KOHNEN/','CR6Series_Metero1min.dat','2018-12-13 12:00','2019-01-22 23:59')
    flux1min['WD_sonic']=np.arctan2(-flux1min.Uy,flux1min.Ux)*180/np.pi
    flux1min['WD']=flux1min.WD_sonic+NoffsetKo
    flux1min['WSpd']=np.sqrt(flux1min.Ux**2+flux1min.Uy**2)
    

    
    fluxcomb_KO=pd.concat([fluxIRGA,fluxTK3_KO],axis=1)
    flux_cleanLE_KO=fluxcomb_KO[['LE','LvE']].dropna()

    # sensible heat
    flux_cleanHS_KO=fluxcomb_KO[['H','HTs']].dropna()


        #comparison daily mean values
    KOHNENdaily=flux_cleanLE_KO.resample('1D').mean()  #LE=IRGASON , LvE=TK3
    KOHNENdailyHs=flux_cleanHS_KO.resample('1D').mean()
    #%% latent heat flux comparison
    fig1ko,(ax1,ax2,ax3)=plt.subplots(3,1)
    ax1.plot(fluxIRGA.LE,color='fuchsia',label='Irgason')
    
    
    ax1.plot(fluxTK3_KO.LvE,color='seagreen',label='TK3')
    ax1.set_title('Comnparison latent heat KOHNEN 18/19')
    ax1.axhline(linestyle=':')
    ax1.set_ylabel(' Wm-2')
    ax1.set_xlim('2018-12-13','2019-01-18')
    
    ax2.plot(flux_cleanLE_KO.LE-flux_cleanLE_KO.LvE,label='IRGA-TK3')
    ax2.set_ylabel(' Wm-2')
    ax2.text(0.2, 0.1, f'mean diff = {np.round(np.mean(flux_cleanLE_KO.LE-flux_cleanLE_KO.LvE),2)}', horizontalalignment='center',verticalalignment='center', transform=ax3.transAxes)
    ax2.text(0.6, 0.1, f'variance diff = {np.round(np.var(flux_cleanLE_KO.LE-flux_cleanLE_KO.LvE),2)}', horizontalalignment='center',verticalalignment='center', transform=ax3.transAxes)
    ax2.set_xlim('2018-12-13','2019-01-18')
    ax2.set_ylim(-5,5)
    ax2.grid(True)
    
    ax3.hist(flux_cleanLE_KO.LE['2018-12-13':'2019-01-18']-flux_cleanLE_KO.LvE['2018-12-13':'2019-01-18'],40,label='diff')
    ax3.set_xlabel(' Wm-2')
    fig1ko.legend(loc='lower right')
    fig1ko.savefig('/Users/swa048/Documents/Data/Kohnen_Data/Kohnen_plots/fluxes/'+'fluxLE_comparison_Kohnen'+'.pdf') 
#%%
    fig2ko,(ax1,ax2,ax3)=plt.subplots(3,1)
    ax1.plot(fluxIRGA.H,color='darkorchid',label='Irgason')
    ax1.set_ylabel('Wm-2')
    ax1.set_xlim('2018-12-13','2019-01-18')
    ax1.set_ylim(-25,35)
    ax2.set_title('13.12.2018-18.01.2019')
    
    #ax2=ax1.twinx()
    ax1.plot(fluxTK3_KO.HTs,color='palegreen',label='TK3')
    ax1.set_title('sensible heat')
    ax1.axhline(linestyle=':')


    
    ax2.plot(flux_cleanHS_KO.H-flux_cleanHS_KO.HTs,label='Irgason-TK3')
    ax2.set_ylabel('Wm-2')    
    ax2.text(0.2, 0.1, f'mean diff = {np.round(np.mean(flux_cleanHS_KO.H-flux_cleanHS_KO.HTs),2)}', horizontalalignment='center',verticalalignment='center', transform=ax3.transAxes)
    ax2.text(0.6, 0.1, f'variance diff = {np.round(np.var(flux_cleanHS_KO.H-flux_cleanHS_KO.HTs),2)}', horizontalalignment='center',verticalalignment='center', transform=ax3.transAxes)
    ax2.set_xlim('2018-12-13','2019-01-18') 
    ax2.set_ylim(-10,10)
    ax3.hist(flux_cleanHS_KO.H['2018-12-13':'2019-01-18']-flux_cleanHS_KO.HTs['2018-12-13':'2019-01-18'],40,label='diff')
    ax3.set_xlabel(' Wm-2')
    fig2ko.legend()
    fig2ko.savefig('/Users/swa048/Documents/Data/Kohnen_Data/Kohnen_plots/fluxes/'+'fluxH_comparison_Kohnen'+'.pdf')  
    
    #%%    
    figdailyko,(ax1,ax2,ax3)=plt.subplots(3,1)
    ax1.plot(KOHNENdaily.LE,color='darkorchid',label='online',MarkerSize=4)
    ax1.set_ylabel('Wm-2')
    ax1.set_xlim('2018-12-13','2019-01-18')
    
    #ax2=ax1.twinx()
    ax1.plot(KOHNENdaily.LvE,color='palegreen',label='TK3',MarkerSize=4)
    ax1.set_title('Latent Heat Daily Means KOHNEN')
    ax1.axhline(linestyle=':')

    
    ax2.plot(KOHNENdaily.LE-KOHNENdaily.LvE,label='online-TK3')
    ax2.set_ylabel(' Wm-2')
    ax2.text(0.2, 0.1, f'mean diff = {np.round(np.mean(KOHNENdaily.LE-KOHNENdaily.LvE),2)}', horizontalalignment='center',verticalalignment='center', transform=ax3.transAxes)
    ax2.text(0.6, 0.1, f'variance diff = {np.round(np.var(KOHNENdaily.LE-KOHNENdaily.LvE),2)}', horizontalalignment='center',verticalalignment='center', transform=ax3.transAxes)
    ax2.set_xlim('2018-12-13','2019-01-18')
    ax2.set_ylim(-1,1)
    ax2.grid(True)
    
    ax3.hist(KOHNENdaily.LE['2018-12-13':'2019-01-18']-KOHNENdaily.LvE['2018-12-13':'2019-01-18'],40,label='diff')
    figdailyko.legend()
    figdailyko.savefig('/Users/swa048/Documents/Data/Kohnen_Data/Kohnen_plots/fluxes/'+'fluxLE_dailycomparison_Kohnen'+'.pdf')   
   
    #%%  sensible heat daily
    figdailykoH,(ax1,ax2,ax3)=plt.subplots(3,1)
    ax1.plot(KOHNENdailyHs.H,color='darkorchid',label='online',MarkerSize=4)
    ax1.set_ylabel('Wm-2')
    ax1.set_xlim('2018-12-13','2019-01-18')
    
    #ax2=ax1.twinx()
    ax1.plot(KOHNENdailyHs.HTs,color='palegreen',label='TK3',MarkerSize=4)
    ax1.set_title('sensible Heat Daily Means KOHNEN')
    ax1.axhline(linestyle=':')

    
    ax2.plot(KOHNENdailyHs.H-KOHNENdailyHs.HTs,label='online-TK3')
    ax2.set_ylabel(' Wm-2')
    ax2.text(0.2, 0.1, f'mean diff = {np.round(np.mean(KOHNENdailyHs.H-KOHNENdailyHs.HTs),2)}', horizontalalignment='center',verticalalignment='center', transform=ax3.transAxes)
    ax2.text(0.6, 0.1, f'variance diff = {np.round(np.var(KOHNENdailyHs.H-KOHNENdailyHs.HTs),2)}', horizontalalignment='center',verticalalignment='center', transform=ax3.transAxes)
    ax2.set_xlim('2018-12-13','2019-01-18')
    ax2.set_ylim(-1,1)
    ax2.grid(True)
    
    ax3.hist(KOHNENdailyHs.H['2018-12-13':'2019-01-18']-KOHNENdailyHs.HTs['2018-12-13':'2019-01-18'],40,label='diff')
    figdailykoH.legend()
    figdailykoH.savefig('/Users/swa048/Documents/Data/Kohnen_Data/Kohnen_plots/fluxes/'+'fluxHs_dailycomparison_Kohnen'+'.pdf')      
#%%    windrose comparison
    windroseKOHNEN=plt.figure(figsize=(24,6))
    
    #Subplot 1: 
    ax1 = windroseKOHNEN.add_subplot(131,projection='windrose')
    ax1.bar(flux1min.WD, flux1min.WSpd, normed=True, opening=1, edgecolor='white',cmap=cm.winter)
    ax1.set_legend()
    ax1.set_title('1min metero data, wind-dir-compass')
    
    #Subplot 2:
    #IRGA had an offset of 60 in program
    ax2 = windroseKOHNEN.add_subplot(132,projection='windrose')
    ax2.bar(fluxIRGA.WD+5, fluxIRGA.WS, normed=True, opening=1, edgecolor='white',cmap=cm.winter)
    ax2.set_legend()
    ax2.set_title('EddyFlux 30min, wind-dir-compass')
    
    #Subplot 3:
    #should already have the offset calculated... might be that 21 ° are missing (delination correction)
    ax3 = windroseKOHNEN.add_subplot(133,projection='windrose')
    ax3.bar(fluxTK3_KO.wind_dir, fluxTK3_KO.u, normed=True, opening=1, edgecolor='white',cmap=cm.winter)
    ax3.set_legend()
    ax3.set_title('TK3 30min, wind-dir-compass')
    windroseKOHNEN.savefig('/Users/swa048/Documents/Data/Kohnen_Data/Kohnen_plots/Wind/'+'windroseKOHNEN'+'.pdf')  

#%% EGRIP 19 28.05.2019-29.07.2019
EGRIP19=1
if EGRIP19:


    
    NoffsetEG19=245 #degrees from North
    
    fluxTK3_cr6=loadfluxperiod('/Users/swa048/forServer/Meteo/EC/2019/CR6/','EGRIP2019_CR6_result_1908161736.csv','2019-05-28 00:00','2019-07-07 15:45')
    #fluxTK3_cr62=loadfluxperiod('/Users/swa048/forServer/Meteo/EC/2019/CR6/','EGRIP2019_CR6_result_1908161736.csv','2019-05-28 00:00','2019-07-07 15:45')
    
    
    #fluxTK3_cr6=pd.concat((fluxTK3_cr6,fluxTK3_cr62),sort=False)#sort=True)
    
    fluxIRGA_EG19=loadCSIfluxperiod('/Users/swa048/forServer/Meteo/EC/2019/CR6/CR6Series_Flux_CSIFormat.dat','2019-05-28 00:00','2019-07-07 16:00')
    
    
    fluxcomb_EG19=pd.concat([fluxIRGA_EG19,fluxTK3_cr6],axis=1)
    removeoutliers(fluxcomb_EG19,'LE',25,-10)
    removeoutliers(fluxcomb_EG19,'LvE',25,-10)
    flux_cleanLE_EG19=fluxcomb_EG19[['LE','LvE']].dropna()

    # sensible heat
    flux_cleanHS_EG19=fluxcomb_EG19[['H','HTs']].dropna()


    #comparison daily mean values
    EGRIP19daily=flux_cleanLE_EG19.resample('1D').mean()  #LE=IRGASON , LvE=TK3
    EGRIP19dailyHs=flux_cleanHS_EG19.resample('1D').mean()
    #%% latent heat flux comparison
    import matplotlib.dates as mdates
    timefmt=mdates.DateFormatter('%y-%m-%d')
    fig1EG19,(ax1,ax2,ax3)=plt.subplots(3,1,figsize=(18,6))
    ax1.plot(fluxIRGA_EG19.LE,color='fuchsia',label='Irgason')

    ax1.plot(fluxTK3_cr6.LvE,color='seagreen',label='TK3')
    ax1.set_title('Comnparison latent heat EGRIP19')
    ax1.axhline(linestyle=':')
    ax1.set_ylabel(' Wm-2')
    ax1.xaxis.set_major_formatter(timefmt)
    ax1.set_ylim(-10,25)
    ax1.legend()
    
    ax2.plot(flux_cleanLE_EG19.LE-flux_cleanLE_EG19.LvE,label='IRGA-TK3')
    ax2.set_ylabel(' Wm-2')
    ax2.text(0.2, 0.1, f'mean diff = {np.round(np.mean(flux_cleanLE_EG19.LE-flux_cleanLE_EG19.LvE),2)}', horizontalalignment='center',verticalalignment='center', transform=ax3.transAxes)
    ax2.text(0.6, 0.1, f'variance diff = {np.round(np.var(flux_cleanLE_EG19.LE-flux_cleanLE_EG19.LvE),2)}', horizontalalignment='center',verticalalignment='center', transform=ax3.transAxes)
    ax2.xaxis.set_major_formatter(timefmt)
    #ax2.set_xlim('2018-12-13','2019-01-18')
    ax2.set_ylim(-5,5)
    ax2.grid(True)
    ax2.legend()
    
    ax3.hist(flux_cleanLE_EG19.LE['2019-05-28':'2019-07-07']-flux_cleanLE_EG19.LvE['2019-05-28':'2019-07-07'],40,label='diff')
    ax3.set_xlabel(' Wm-2')
    ax3.legend() #loc='lower right'
    fig1EG19.savefig('/Users/swa048/Documents/OneDrive/Dokumente/EGRIP/EGRIP2019/EGRIP_plots2019/flux/'+'fluxLE_comparison_EGRIP19'+'.pdf') 