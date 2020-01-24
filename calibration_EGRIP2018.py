# -*- coding: utf-8 -*-
"""
Created on Wed Sep 26 10:37:48 2018

calibration TC and wind of EC system

@author: Sonja
"""

import numpy as np
import matplotlib.pyplot as plt
import os
import matplotlib.dates as mdates
#from pandas.compat import StringIO
from io import BytesIO
import pandas as pd
import datetime as dt
import numpy.polynomial.polynomial as poly
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

from loadandconvert import loadperiod,exportdf


#%% before season
#calibration period 10.05.18 11:40 - 12.05.18 7:40


ECfile_start='/Users/swa048/forServer/Meteo/Ambient/2018/EGRIP/calibration/beginningofseason/CR3000_ec_scfd_cali.dat'
ATMfile_start='/Users/swa048/forServer/Meteo/Ambient/2018/EGRIP/calibration/beginningofseason/CR3000_Table_ae_cali.dat'

path_plot='/Users/swa048/forServer/Meteo/Ambient/2018/EGRIP/calibration/'
plt.close('all')
plotting=1

#%%
#and import 10 min averages as dataframe
ECs=pd.read_csv(ECfile_start,index_col=0,skiprows=[0,2,3],parse_dates=True,na_values=['NAN'] )
ECs.rename(columns={'RECORD': 'DOY','Ts_Avg':'Ts'}, inplace=True)
t0=np.datetime64(str(ECs.index[0])[0:4]+'-01-01')   #used to be periodstart[0:4]
ECs.loc[:,'DOY'] = pd.Series((ECs.index-t0)/dt.timedelta(days=1), index=ECs.index)    
#df.rename(columns={'ln_vh_Avg': 'hum'}, inplace=True)
#df.loc[:,'hum']=pd.Series(df.hum/(-0.105),index=df.index)

namesall=['timestamp','DOY', 'TC1', 'TC2', 'TC3', 'TC4', 'TC5', 'TC6', 'TC7', 'WS_1','WS_2', 'WS_3', 'WS_1max', 'WS_2max','WS_3max', 'BattV', 'PTemp_C']

#import atmospheric parameters 
ATMs=pd.read_csv(ATMfile_start,index_col=0,skiprows=[0,1,2,3],parse_dates=True,infer_datetime_format=True,na_values=['NAN'],names=namesall )
t0=np.datetime64(str(ATMs.index[0])[0:4]+'-01-01')   #used to be periodstart[0:4]
ATMs.loc[:,'DOY'] = pd.Series((ATMs.index-t0)/dt.timedelta(days=1), index=ATMs.index) 
ATMs_cali=ATMs.copy()
##sozusagen als template
#ATMs_cali=pd.read_csv(path_start+ATMfile_start,index_col=0,skiprows=[0,2,3],parse_dates=True,na_values=['NAN'],dtype={'RECORD': np.int64, 'FWTC_C_Avg(1)': np.float64, 'FWTC_C_Avg(2)': np.float64, 'FWTC_C_Avg(3)': np.float64, 'FWTC_C_Avg(4)': np.float64, 'FWTC_C_Avg(5)': np.float64, 'FWTC_C_Avg(6)': np.float64, 'FWTC_C_Avg(7)': np.float64, 'WS_ms_Avg(1)': np.float64, 'WS_ms_Avg(2)': np.float64, 'WS_ms_Avg(3)': np.float64, 'WS_ms_Max(1)': np.float64, 'WS_ms_Max(2)': np.float64, 'WS_ms_Max(3)': np.float64, 'BattV': np.float64, 'PTemp_C': np.float64},usecols=['TIMESTAMP','RECORD', 'FWTC_C_Avg(1)', 'FWTC_C_Avg(2)', 'FWTC_C_Avg(3)',  'FWTC_C_Avg(5)', 'FWTC_C_Avg(6)', 'FWTC_C_Avg(7)', 'WS_ms_Avg(1)', 'WS_ms_Avg(2)', 'WS_ms_Avg(3)'] )
#ATMs_cali.rename(columns={'RECORD': 'DOY','FWTC_C_Avg(1)':'TC1','FWTC_C_Avg(2)':'TC2','FWTC_C_Avg(3)':'TC3','FWTC_C_Avg(4)':'TC4','FWTC_C_Avg(5)':'TC5','FWTC_C_Avg(6)':'TC6','FWTC_C_Avg(7)':'TC7',}, inplace=True)
#ATMs_cali.loc[:,'DOY'] = pd.Series(ATMs_cali.index.dayofyear, index=ATMs_cali.index)
# =============================================================================
# good way of using read_csv
# df = pd.read_csv(StringIO(csv),
#         header=0,
#         index_col=["date", "loc"], 
#         usecols=["date", "loc", "x"],
#         parse_dates=["date"])
# =============================================================================





#%%
if plotting:
    plt.figure(figsize=(16,8))
    plt.plot(ATMs.TC5,ECs.Ts,'o',color='blueviolet',markersize=5) #2cm
    plt.plot(ATMs.TC6,ECs.Ts,'o',color='crimson',markersize=5) #5cm
    plt.plot(ATMs.TC7,ECs.Ts,'o',color='indianred',markersize=5) #10cm
    plt.plot(ATMs.TC1,ECs.Ts,'o',color='goldenrod',markersize=5) #30cm
    plt.plot(ATMs.TC2,ECs.Ts,'o',color='yellowgreen',markersize=5) #60cm
    plt.plot(ATMs.TC3,ECs.Ts,'o',color='darkcyan',markersize=5) #180cm
    #something seems to be have happen with TC4
    plt.xlabel('sonic temperature C')
    plt.ylabel('TC temperature C')
    plt.show()


#calibrate all Tc against Ts and put them in ATMs_cali dataframe
coeffsTC=pd.DataFrame(index=ATMs.columns[[1,2,3,5,6,7]],columns=['c','m'])

for tc in ATMs.columns[[1,2,3,5,6,7]]:
    coefs,rr  = poly.polyfit(ATMs[tc], ECs.Ts, 1,full=True)  #coeff, r  with r= [sumofresiduals, rank, singular_values, rcond] : list
    coeffsTC.loc[tc,:]=coefs
    ATMs_cali[tc] = poly.polyval(ATMs[tc], coefs)
    if plotting:
        if tc == ATMs.columns[5]:
            fig,ax1 = plt.subplots(figsize=(16,8))
            ax1.plot(ATMs[tc],ECs.Ts,'o',color='goldenrod',markersize=5) #30cm
            ax1.plot(ATMs[tc], ATMs_cali[tc],color='navy')
            ax1.set_xlabel('TC temp C')
            ax1.set_ylabel('T sonic temp C')
            ax1.set_title(tc+'start')
            ax1.legend(['data','fit'])
            
            ax2=fig.add_axes([0.6,0.15,0.25,0.25],title='Residuals      sum= '+str(rr[0]))
            ax2.hist(ECs.Ts-ATMs_cali[tc])
            ax2.set_xlabel=('count')
            plt.figtext(0.2,0.65,('coefs:  m='+str(coefs[1])+' c='+str(coefs[0])))
            fig.set_tight_layout
            plt.close(fig)


#%% after season

#calibration period 03.08.18-07.8.18


ECfile_end='/Users/swa048/forServer/Meteo/Ambient/2018/EGRIP/calibration/endofseason/CR3000_ec_scfd_calibration_end.dat'
ATMfile_end='/Users/swa048/forServer/Meteo/Ambient/2018/EGRIP/calibration/endofseason/CR3000_Table_ae_calibration_end.dat'



namesall=['timestamp','DOY', 'TC1', 'TC2', 'TC3', 'TC4', 'TC5', 'TC6', 'TC7', 'WS_1','WS_2', 'WS_3', 'WS_1max', 'WS_2max','WS_3max', 'BattV', 'PTemp_C']

#import atmospheric parameters 
ATMe=pd.read_csv(ATMfile_end,index_col=0,skiprows=[0,1,2,3],parse_dates=True,infer_datetime_format=True,na_values=['NAN'],names=namesall )
t0=np.datetime64(str(ATMe.index[0])[0:4]+'-01-01')   #used to be periodstart[0:4]
ATMe.loc[:,'DOY'] = pd.Series((ATMe.index-t0)/dt.timedelta(days=1), index=ATMe.index) 
#auch wieder als vorlge für die calibrierten daten
ATMe_cali=ATMe.copy()


#and import 10 min averages as dataframe
ECeall=pd.read_csv(ECfile_end,index_col=0,parse_dates=True,na_values=['NAN'] )
ECeall.rename(columns={'RECORD': 'DOY','Ts_Avg':'Ts'}, inplace=True)
t0=np.datetime64(str(ECeall.index[0])[0:4]+'-01-01')   #used to be periodstart[0:4]
ECeall.loc[:,'DOY'] = pd.Series((ECeall.index-t0)/dt.timedelta(days=1), index=ECeall.index)    
#df.rename(columns={'ln_vh_Avg': 'hum'}, inplace=True)
#df.loc[:,'hum']=pd.Series(df.hum/(-0.105),index=df.index)

ECe=ECeall.loc[ATMe.index[0]:ATMe.index[-1],:]


#%%
#calibrate all Tc against Ts and put them in ATMs_cali dataframe
#TC4 and TC5 are broken and do not need to be 
coeffeTC=pd.DataFrame(index=ATMe.columns[[1,2,3,6,7]],columns=['c','m'])

for tce,tcs,ixs in zip(ATMe.columns[[1,2,3,6,7]],ATMs.columns[[1,2,3,6,7]],[0,1,2,4,5]):  #the first two idex with "names" and the last (ixs) with indeces
    coefe,rre  = poly.polyfit(ATMe[tce], ECe.Ts, 1,full=True)  #coeff, r  with r= [sumofresiduals, rank, singular_values, rcond] : list
    coeffeTC.loc[tce,:]=coefe
    ATMe_cali[tce] = poly.polyval(ATMe[tce], coefe)
    
    
if plotting:
    for tce in coeffeTC.index:
        fig,ax = plt.subplots(1,1,figsize=(8,12))
        ax.plot(ATMe[tce],ECe.Ts,'o',color='goldenrod',markersize=5) #30cm
        ax.plot(ATMe[tce], ATMe_cali[tce],color='navy')
        ax.plot(ATMs[tce],ECs.Ts,'o',color='cyan',markersize=5) #30cm
        ax.plot(ATMs[tce], ATMs_cali[tce],color='yellowgreen')
        ax.set_xlabel('TC temp C')
        ax.set_ylabel('T sonic temp C')
        ax.set_title(tce+'end')
        ax.legend(['data','fit'])
        ax.text(-20,-10,('coefs:  m='+str(coeffeTC.loc[tce,'m'])+' c='+str(coeffeTC.loc[tce,'c'])))
        
        # plot in subplot 1
        ax02= inset_axes(ax,height="25%", # set height
                              width="25%", # and width
                              loc=4)
        ax02.hist(ECe.Ts-ATMe_cali[tce])
        ax02.set_xlabel=('count')
        ax02.set_title('Residuals  sum= '+str(rr[0]))
        
    
        
        #subplot 2

#        ax.set_xlabel('TC temp C')
#        ax.set_ylabel('T sonic temp C')
#        ax.set_title(tcs+'start')
#        ax.legend(['data','fit'])
#        ax.text(-32,-22,('coefs:  m='+str(coeffsTC.loc[tcs,'m'])+' c='+str(coeffsTC.loc[tcs,'c'])))
#        
#        #plot in subplot 2
#        ax12=inset_axes(ax,height="25%", # set height
#                              width="25%", # and width
#                              loc=4)
#        ax12.hist(ECs.Ts-ATMs_cali[tcs])
#        ax12.set_title='Residuals  sum= '+str(rr[0])
#        ax12.set_xlabel=('count')
        fig.set_tight_layout

        fig.savefig(path_plot+str(tce)+'_valuesTCcali'+'.png')
        plt.close(fig)

#%% combined calibration periods

ATMcombi=pd.concat([ATMs,ATMe],axis=0)  #TC4 & 5 are broken; not usable 
ECcombi=pd.concat([ECs,ECe],axis=0)

coeffTCcombi=pd.DataFrame(index=ATMe.columns[[1,2,3,6,7]],columns=['c','m'])


for tc in ATMcombi.columns[[1,2,3,6,7]]:  #counter,value
    coeffTCcombi.loc[tc,:] = poly.polyfit(ATMcombi[tc],ECcombi.Ts,  1,full=False) #polyfit(x,y,degree)

#add the calibration values where we only have start calibration values    
# TC4 is broken from early on
coeffTCcombi.loc['TC5',:]=coeffsTC.loc['TC5',:]
#%%
#to check if calibration is visually okay
TCcalibrated_combi=pd.DataFrame(data=np.zeros_like(ATMcombi[ATMcombi.columns[[1,2,3,6,7]]]),index=ATMcombi.index,columns=ATMe.columns[[1,2,3,6,7]])  # 1st row as the column names

#load full season 10 min data
path_fullseason='/Users/swa048/forServer/Meteo/Ambient/2018/EGRIP/data/CR3000_Table_ae.dat'
ATM_fullseason=pd.read_csv(path_fullseason,index_col=0,skiprows=[0,1,2,3],parse_dates=True,infer_datetime_format=True,na_values=['NAN'],names=namesall )
t0=np.datetime64(str(ATM_fullseason.index[0])[0:4]+'-01-01')   #used to be periodstart[0:4]
ATM_fullseason.loc[:,'DOY'] = pd.Series((ATM_fullseason.index-t0)/dt.timedelta(days=1), index=ATM_fullseason.index) 

#load full season 1s data for calibration
path_1s='/Users/swa048/forServer/Meteo/Ambient/2018/EGRIP/data/CR3000_Table_wT1s.dat'
ATM_1s=pd.read_csv(path_1s,index_col=0,skiprows=[0,1,2,3],parse_dates=True,infer_datetime_format=True,na_values=['NAN'],names=namesall )
t0=np.datetime64(str(ATM_1s.index[0])[0:4]+'-01-01')   #used to be periodstart[0:4]
ATM_1s.loc[:,'DOY'] = pd.Series((ATM_1s.index-t0)/dt.timedelta(days=1), index=ATM_1s.index) 


#load full season 10 min data from Irgason
ECall=pd.read_csv('/Users/swa048/forServer/Meteo/EC/2018/EGRIP/onlinefluxes10min/CR3000_ec_scfd.dat',index_col=0,skiprows=[0,2,3],parse_dates=True,na_values=['NAN'] )
ECall.rename(columns={'RECORD': 'DOY','Ts_Avg':'Ts'}, inplace=True)
t0=np.datetime64(str(ECall.index[0])[0:4]+'-01-01')   #used to be periodstart[0:4]
ECall.loc[:,'DOY'] = pd.Series((ECall.index-t0)/dt.timedelta(days=1), index=ECall.index) 
#%%  EIGETNLICHE CALIBRATION

#apply calibration to full season TC record which is named ATM_fullseason and save it to dataframe ATMcalibrated_all which can be exported
ATMcalibrated_all=ATM_fullseason.copy()  # 1st row as the column names
ATM_1scali=ATM_1s.copy()  # 1st row as the column names

for tcall in coeffTCcombi.index:  #counter,value     TC 5 is the beginning calibration values
    TCcalibrated_combi.loc[:,tcall] = poly.polyval(ATMcombi.loc[:,tcall], coeffTCcombi.loc[tcall,:])  #eigetnliche calibrierung
    ATMcalibrated_all.loc[:,tcall] = poly.polyval(ATM_fullseason.loc[:,tcall], coeffTCcombi.loc[tcall,:])  #eigetnliche calibrierung
    ATM_1scali.loc[:,tcall] = poly.polyval(ATM_1s.loc[:,tcall], coeffTCcombi.loc[tcall,:])  #eigetnliche calibrierung
#%%    
if plotting:
   for tcall in coeffTCcombi.index: 
        fig,ax1 = plt.subplots(figsize=(16,8))
        ax2=fig.add_axes([0.6,0.15,0.25,0.25],title='Residuals      sum= '+str(rr[0]))
        if tcall == 'TC5':
            ax1.plot(ATMs[tcall],ECs.Ts,'o',color='goldenrod',markersize=5) #30cm
            ax1.plot(ATMs[tcall], ATMs_cali[tcall],color='navy')
            ax1.set_title(tcall+' start')
            ax2.hist(ECs.Ts-ATMs_cali[tcall])
            plt.figtext(0.2,0.65,('coefs:  m='+str(coeffsTC.loc[tcall,'m'])+' c='+str(coeffsTC.loc[tcall,'c'])))
        else:    
            ax1.plot(ATMcombi.loc[:,tcall],ECcombi.Ts,'o',color='goldenrod',markersize=5) 
            ax1.plot(ATMcombi.loc[:,tcall], TCcalibrated_combi.loc[:,tcall],color='navy')
            ax1.set_title(tcall+' combi')
            ax2.hist((ECcombi.Ts-TCcalibrated_combi.loc[:,tcall]).dropna())
            plt.figtext(0.2,0.65,('coefs:  m='+str(coeffTCcombi.loc[tcall,'m'])+' c='+str(coeffTCcombi.loc[tcall,'c'])))
        
        
        ax1.set_xlabel('TC temp C')
        ax1.set_ylabel('T sonic temp C')
        ax1.legend(['data','fit'])
        ax2.set_xlabel=('count')
        
        fig.set_tight_layout
        fig.savefig(path_plot+str(tcall)+'_combinedTCcali'+'.png')
        plt.close(fig)
#%%
#calibrate wind sensors

#middle wind sensor should not be calibrated since it is out of order --> WS_2
#only have WS and wind direction in one calibration dataframe

wind_combi=ATMcombi.loc[:,['WS_1', 'WS_2','WS_3']] #WS_1 = 30cm; WS_3=180cm
windEC_combi=ECcombi.loc[:,['Ux_Avg','Uy_Avg', 'Uz_Avg','horiz_wind_spd','result_wind_spd']]

windcalibrated_combi=pd.DataFrame(data=np.zeros_like(wind_combi),index=wind_combi.index,columns=wind_combi.columns) 
#%%
figwind,((ax1,ax2),(ax3,ax4),(ax5,ax6))=plt.subplots(3,2,figsize=(16,16))
ax1.plot(wind_combi.loc[ATMs.index[0]:ATMs.index[-1],'WS_1'],windEC_combi.horiz_wind_spd[ATMs.index[0]:ATMs.index[-1]],'*',label='start')
ax1.plot(wind_combi.loc[ATMe.index[0]:ATMe.index[-1],'WS_1'],windEC_combi.horiz_wind_spd[ATMe.index[0]:ATMe.index[-1]],'*',color='orange',label='end')
ax2.plot(wind_combi.loc[ATMs.index[0]:ATMs.index[-1],'WS_1'],windEC_combi.result_wind_spd[ATMs.index[0]:ATMs.index[-1]],'*',label='start')
ax2.plot(wind_combi.loc[ATMe.index[0]:ATMe.index[-1],'WS_1'],windEC_combi.result_wind_spd[ATMe.index[0]:ATMe.index[-1]],'*',color='orange',label='end')

ax1.set_ylabel('sonic anemometer')


###
ax3.plot(wind_combi.loc[ATMs.index[0]:ATMs.index[-1],'WS_2'],windEC_combi.horiz_wind_spd[ATMs.index[0]:ATMs.index[-1]],'*',label='start')
ax3.plot(wind_combi.loc[ATMe.index[0]:ATMe.index[-1],'WS_2'],windEC_combi.horiz_wind_spd[ATMe.index[0]:ATMe.index[-1]],'*',color='orange',label='end')
ax4.plot(wind_combi.loc[ATMs.index[0]:ATMs.index[-1],'WS_2'],windEC_combi.result_wind_spd[ATMs.index[0]:ATMs.index[-1]],'*',label='start')
ax4.plot(wind_combi.loc[ATMe.index[0]:ATMe.index[-1],'WS_2'],windEC_combi.result_wind_spd[ATMe.index[0]:ATMe.index[-1]],'*',color='orange',label='end')

ax3.set_ylabel('sonic anemometer')


###
ax5.plot(wind_combi.loc[ATMs.index[0]:ATMs.index[-1],'WS_3'],windEC_combi.horiz_wind_spd[ATMs.index[0]:ATMs.index[-1]],'*',label='start')
ax5.plot(wind_combi.loc[ATMe.index[0]:ATMe.index[-1],'WS_3'],windEC_combi.horiz_wind_spd[ATMe.index[0]:ATMe.index[-1]],'*',color='orange',label='end')
ax6.plot(wind_combi.loc[ATMs.index[0]:ATMs.index[-1],'WS_3'],windEC_combi.result_wind_spd[ATMs.index[0]:ATMs.index[-1]],'*',label='start')
ax6.plot(wind_combi.loc[ATMe.index[0]:ATMe.index[-1],'WS_3'],windEC_combi.result_wind_spd[ATMe.index[0]:ATMe.index[-1]],'*',color='orange',label='end')
ax5.set_xlabel('cup anemometer')
ax5.set_ylabel('sonic anemometer')
ax6.set_xlabel('cup anemometer')


ax1.set_title('WS1 horiz_wind_spd')
ax2.set_title('WS1 result_wind_spd')

ax3.set_title('WS2 horiz_wind_spd - broken')
ax4.set_title('WS2 result_wind_spd - broken')

ax5.set_title('WS3 horiz_wind_spd')
ax6.set_title('WS3 result_wind_spd')
ax1.legend()
ax2.legend()
figwind.tight_layout()
figwind.savefig(path_plot+'_valueswindcali'+'.png')
#%%
plt.figure()
plt.plot(windEC_combi.result_wind_spd-windEC_combi.horiz_wind_spd)

# --> just go with the horizontal wind speed, there is no big difference anyways
#%% calibration of the wind_1 which is the one at 30 cm


coeff_wind=pd.DataFrame(index=['WS_1','WS_3'],columns=['c','m'])


for wind in coeff_wind.index:  #counter,value
    coeff_wind.loc[wind,:] = poly.polyfit(wind_combi.loc[:,wind],windEC_combi.horiz_wind_spd,  1,full=False) #polyfit(x,y,degree)

for windcall in coeff_wind.index:  #counter,value
    windcalibrated_combi.loc[:,windcall] = poly.polyval(wind_combi.loc[:,windcall], coeff_wind.loc[windcall,:])  #eigetnliche calibrierung
    ATMcalibrated_all.loc[:,windcall] = poly.polyval(ATM_fullseason.loc[:,windcall], coeff_wind.loc[windcall,:])  #eigetnliche calibrierung


if plotting:
   for windcall in coeff_wind.index: 
        fig,ax1 = plt.subplots(figsize=(16,8))
        ax1.plot(wind_combi.loc[:,windcall],windEC_combi.horiz_wind_spd,'o',color='darkgreen',markersize=5) 
        ax1.plot(wind_combi.loc[:,windcall], windcalibrated_combi.loc[:,windcall],color='crimson')
        ax1.set_xlabel('cup anemometer')
        ax1.set_ylabel('Sonic / calibrated')
        ax1.set_title(windcall+' combi')
        ax1.legend(['data','fit'])
            
        ax2=fig.add_axes([0.6,0.15,0.25,0.25],title='Residuals      sum ')
        ax2.hist(windEC_combi.horiz_wind_spd-windcalibrated_combi.loc[:,windcall])
        ax2.set_xlabel=('count')
        plt.figtext(0.2,0.65,('coefs:  m='+str(coeff_wind.loc[windcall,'m'])+' c='+str(coeff_wind.loc[windcall,'c'])))
        fig.set_tight_layout
        fig.savefig(path_plot+str(windcall)+'_windcali'+'.png')
        plt.close(fig)
        
#%%
plt.figure()
plt.plot(ATMcalibrated_all.loc[:,'WS_1'],label='cup 30cm')       
plt.plot(ATMcalibrated_all.loc[:,'WS_3'],label='cup 180cm')   
plt.plot(ECall.loc[:,'horiz_wind_spd'],label='CSAT')
plt.legend()  
        
#%%
#export calibrated 10 min data

ATMcalibrated_all.drop(columns=['TC4','WS_2','WS_1max', 'WS_2max', 'WS_3max', 'BattV', 'PTemp_C'],inplace=True)
ATMcalibrated_all.rename(inplace=True,columns={'TC1': 'TC_30cm',
                                                'TC2': 'TC_80cm',
                                                'TC3': 'TC_180cm',
                                                'TC5': 'TC_2cm',
                                                'TC6': 'TC_5cm',
                                                'TC7': 'TC_10cm',
                                                'WS_1':'WS_30cm',
                                                'WS_3':'WS_180cm'})
ATM_1scali.drop(columns=['TC4','WS_2','WS_1max', 'WS_2max', 'WS_3max', 'BattV', 'PTemp_C'],inplace=True)
ATM_1scali.rename(inplace=True,columns={'TC1': 'TC_30cm',
                                                'TC2': 'TC_80cm',
                                                'TC3': 'TC_180cm',
                                                'TC5': 'TC_2cm',
                                                'TC6': 'TC_5cm',
                                                'TC7': 'TC_10cm',
                                                'WS_1':'WS_30cm',
                                                'WS_3':'WS_180cm'})

exportdf(ATMcalibrated_all,'/Users/swa048/forServer/Meteo/Ambient/2018/EGRIP/data/','MeteoEGRIP18_10min_calibrated')
exportdf(ATM_1scali,'/Users/swa048/forServer/Meteo/Ambient/2018/EGRIP/data/','MeteoEGRIP18_1s_calibrated')
#%%
#check if import works
dfdummy=pd.read_csv('/Users/swa048/forServer/Meteo/Ambient/2018/EGRIP/data/MeteoEGRIP18_10min_calibrated.txt',index_col=0,parse_dates=True,na_values=['NAN'])