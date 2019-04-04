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
from pandas.compat import StringIO
from io import BytesIO
import pandas as pd
from datetime import datetime
import numpy.polynomial.polynomial as poly
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

from loadandconvert import loadperiod,exportdf


#%% before season
#calibration period 10.05.18 11:40 - 12.05.18 7:40

path_start='/Volumes/U/EC_DATA/Calibration_beginning_season/'
ECfile_start='CR3000_ec_scfd_cali.dat'
ATMfile_start='CR3000_Table_ae_cali.dat'

path_plot='/Volumes/U/EC_processed/plots/daily/'
plt.close('all')
plotting=0


#and import 10 min averages as dataframe
ECs=pd.read_csv(path_start+ECfile_start,index_col=0,skiprows=[0,2,3],parse_dates=True,na_values=['NAN'],dtype={'RECORD': np.int64, 'LE': np.float64, 'Hs': np.float64, 'tau': np.float64, 'u_star': np.float64, 'Ux_Ux': np.float64, 'Ux_Uy': np.float64, 'Ux_Uz': np.float64, 'Ux_ln_vh': np.float64, 'Ux_Ts': np.float64, 'Uy_Uy': np.float64, 'Uy_Uz': np.float64, 'Uy_ln_vh': np.float64, 'Uy_Ts': np.float64, 'Uz_Uz': np.float64, 'Uz_ln_vh': np.float64, 'Uz_Ts': np.float64, 'ln_vh_ln_vh': np.float64, 'ln_vh_Ts': np.float64, 'Ts_Ts': np.float64, 'Ux_Avg': np.float64, 'Uy_Avg': np.float64, 'Uz_Avg': np.float64, 'ln_vh_Avg': np.float64, 'Ts_Avg': np.float64, 'horiz_wind_spd': np.float64, 'result_wind_spd': np.float64, 'wind_dir': np.float64, 'wind_dir_std_dev': np.float64, 'wnd_dir_compass': np.float64, 'vh_Avg': np.float64, 'n_Tot': np.float64, 'del_T_f_Tot': np.float64, 'track_f_Tot': np.float64, 'amp_h_f_Tot': np.float64, 'amp_l_f_Tot': np.float64} )
ECs.rename(columns={'RECORD': 'DOY','Ts_Avg':'Ts'}, inplace=True)
ECs.loc[:,'DOY'] = pd.Series(ECs.index.dayofyear, index=ECs.index)
#df.rename(columns={'ln_vh_Avg': 'hum'}, inplace=True)
#df.loc[:,'hum']=pd.Series(df.hum/(-0.105),index=df.index)


#import atmospheric parameters 
ATMs=pd.read_csv(path_start+ATMfile_start,index_col=0,skiprows=[0,2,3],parse_dates=True,na_values=['NAN'],dtype={'RECORD': np.int64, 'FWTC_C_Avg(1)': np.float64, 'FWTC_C_Avg(2)': np.float64, 'FWTC_C_Avg(3)': np.float64, 'FWTC_C_Avg(4)': np.float64, 'FWTC_C_Avg(5)': np.float64, 'FWTC_C_Avg(6)': np.float64, 'FWTC_C_Avg(7)': np.float64, 'WS_ms_Avg(1)': np.float64, 'WS_ms_Avg(2)': np.float64, 'WS_ms_Avg(3)': np.float64, 'WS_ms_Max(1)': np.float64, 'WS_ms_Max(2)': np.float64, 'WS_ms_Max(3)': np.float64, 'BattV': np.float64, 'PTemp_C': np.float64} )
ATMs.rename(columns={'RECORD': 'DOY','FWTC_C_Avg(1)':'TC1','FWTC_C_Avg(2)':'TC2','FWTC_C_Avg(3)':'TC3','FWTC_C_Avg(4)':'TC4','FWTC_C_Avg(5)':'TC5','FWTC_C_Avg(6)':'TC6','FWTC_C_Avg(7)':'TC7',}, inplace=True)
ATMs.loc[:,'DOY'] = pd.Series(ATMs.index.dayofyear, index=ATMs.index)

#sozusagen als template
ATMs_cali=pd.read_csv(path_start+ATMfile_start,index_col=0,skiprows=[0,2,3],parse_dates=True,na_values=['NAN'],dtype={'RECORD': np.int64, 'FWTC_C_Avg(1)': np.float64, 'FWTC_C_Avg(2)': np.float64, 'FWTC_C_Avg(3)': np.float64, 'FWTC_C_Avg(4)': np.float64, 'FWTC_C_Avg(5)': np.float64, 'FWTC_C_Avg(6)': np.float64, 'FWTC_C_Avg(7)': np.float64, 'WS_ms_Avg(1)': np.float64, 'WS_ms_Avg(2)': np.float64, 'WS_ms_Avg(3)': np.float64, 'WS_ms_Max(1)': np.float64, 'WS_ms_Max(2)': np.float64, 'WS_ms_Max(3)': np.float64, 'BattV': np.float64, 'PTemp_C': np.float64},usecols=['TIMESTAMP','RECORD', 'FWTC_C_Avg(1)', 'FWTC_C_Avg(2)', 'FWTC_C_Avg(3)',  'FWTC_C_Avg(5)', 'FWTC_C_Avg(6)', 'FWTC_C_Avg(7)', 'WS_ms_Avg(1)', 'WS_ms_Avg(2)', 'WS_ms_Avg(3)'] )
ATMs_cali.rename(columns={'RECORD': 'DOY','FWTC_C_Avg(1)':'TC1','FWTC_C_Avg(2)':'TC2','FWTC_C_Avg(3)':'TC3','FWTC_C_Avg(4)':'TC4','FWTC_C_Avg(5)':'TC5','FWTC_C_Avg(6)':'TC6','FWTC_C_Avg(7)':'TC7',}, inplace=True)
ATMs_cali.loc[:,'DOY'] = pd.Series(ATMs_cali.index.dayofyear, index=ATMs_cali.index)
# =============================================================================
# good way of using read_csv
# df = pd.read_csv(StringIO(csv),
#         header=0,
#         index_col=["date", "loc"], 
#         usecols=["date", "loc", "x"],
#         parse_dates=["date"])
# =============================================================================


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
coeffsTC=[]

for tc in ATMs.columns[[1,2,3,5,6,7]]:
    coefs,rr  = poly.polyfit(ATMs[tc], ECs.Ts, 1,full=True)  #coeff, r  with r= [sumofresiduals, rank, singular_values, rcond] : list
    coeffsTC.append(np.append(tc+'start',coefs).reshape((1,3)))
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
            fig.show()

#%% after season

#calibration period 03.08.18-07.8.18

path_end='/Volumes/U/EC_DATA/20180807_endofcalibration/'
ECfile_end='CR3000_ec_scfd_calibration_end.dat'
ATMfile_end='CR3000_Table_ae_calibration_end.dat'

path_plot='/Volumes/U/EC_processed/plots/daily/'



#and import 10 min averages as dataframe
ECe=pd.read_csv(path_end+ECfile_end,index_col=0,parse_dates=True,na_values=['NAN',-7999],dtype={'RECORD': np.int64, 'LE': np.float64, 'Hs': np.float64, 'tau': np.float64, 'u_star': np.float64, 'Ux_Ux': np.float64, 'Ux_Uy': np.float64, 'Ux_Uz': np.float64, 'Ux_ln_vh': np.float64, 'Ux_Ts': np.float64, 'Uy_Uy': np.float64, 'Uy_Uz': np.float64, 'Uy_ln_vh': np.float64, 'Uy_Ts': np.float64, 'Uz_Uz': np.float64, 'Uz_ln_vh': np.float64, 'Uz_Ts': np.float64, 'ln_vh_ln_vh': np.float64, 'ln_vh_Ts': np.float64, 'Ts_Ts': np.float64, 'Ux_Avg': np.float64, 'Uy_Avg': np.float64, 'Uz_Avg': np.float64, 'ln_vh_Avg': np.float64, 'Ts_Avg': np.float64, 'horiz_wind_spd': np.float64, 'result_wind_spd': np.float64, 'wind_dir': np.float64, 'wind_dir_std_dev': np.float64, 'wnd_dir_compass': np.float64, 'vh_Avg': np.float64, 'n_Tot': np.float64, 'del_T_f_Tot': np.float64, 'track_f_Tot': np.float64, 'amp_h_f_Tot': np.float64, 'amp_l_f_Tot': np.float64} )
ECe.rename(columns={'RECORD': 'DOY','Ts_Avg':'Ts'}, inplace=True)
ECe.loc[:,'DOY'] = pd.Series(ECe.index.dayofyear, index=ECe.index)
#df.rename(columns={'ln_vh_Avg': 'hum'}, inplace=True)
#df.loc[:,'hum']=pd.Series(df.hum/(-0.105),index=df.index)


#import atmospheric parameters 
ATMe=pd.read_csv(path_end+ATMfile_end,index_col=0,parse_dates=True,na_values=['NAN',-7999],dtype={'RECORD': np.int64, 'FWTC_C_Avg(1)': np.float64, 'FWTC_C_Avg(2)': np.float64, 'FWTC_C_Avg(3)': np.float64, 'FWTC_C_Avg(4)': np.float64, 'FWTC_C_Avg(5)': np.float64, 'FWTC_C_Avg(6)': np.float64, 'FWTC_C_Avg(7)': np.float64, 'WS_ms_Avg(1)': np.float64, 'WS_ms_Avg(2)': np.float64, 'WS_ms_Avg(3)': np.float64, 'WS_ms_Max(1)': np.float64, 'WS_ms_Max(2)': np.float64, 'WS_ms_Max(3)': np.float64, 'BattV': np.float64, 'PTemp_C': np.float64} )
ATMe.rename(columns={'RECORD': 'DOY','FWTC_C_Avg(1)':'TC1','FWTC_C_Avg(2)':'TC2','FWTC_C_Avg(3)':'TC3','FWTC_C_Avg(4)':'TC4','FWTC_C_Avg(5)':'TC5','FWTC_C_Avg(6)':'TC6','FWTC_C_Avg(7)':'TC7',}, inplace=True)
ATMe.loc[:,'DOY'] = pd.Series(ATMe.index.dayofyear, index=ATMe.index)

#auch wieder als vorlge für die calibrierten daten
ATMe_cali=pd.read_csv(path_end+ATMfile_end,index_col=0,parse_dates=True,na_values=['NAN',-7999],dtype={'RECORD': np.int64, 'FWTC_C_Avg(1)': np.float64, 'FWTC_C_Avg(2)': np.float64, 'FWTC_C_Avg(3)': np.float64, 'FWTC_C_Avg(4)': np.float64, 'FWTC_C_Avg(5)': np.float64, 'FWTC_C_Avg(6)': np.float64, 'FWTC_C_Avg(7)': np.float64, 'WS_ms_Avg(1)': np.float64, 'WS_ms_Avg(2)': np.float64, 'WS_ms_Avg(3)': np.float64, 'WS_ms_Max(1)': np.float64, 'WS_ms_Max(2)': np.float64, 'WS_ms_Max(3)': np.float64, 'BattV': np.float64, 'PTemp_C': np.float64},usecols=['TIMESTAMP','RECORD', 'FWTC_C_Avg(1)', 'FWTC_C_Avg(2)', 'FWTC_C_Avg(3)', 'FWTC_C_Avg(5)','FWTC_C_Avg(6)', 'FWTC_C_Avg(7)', 'WS_ms_Avg(1)', 'WS_ms_Avg(2)', 'WS_ms_Avg(3)'] )
ATMe_cali.rename(columns={'RECORD': 'DOY','FWTC_C_Avg(1)':'TC1','FWTC_C_Avg(2)':'TC2','FWTC_C_Avg(3)':'TC3','FWTC_C_Avg(4)':'TC4','FWTC_C_Avg(5)':'TC5','FWTC_C_Avg(6)':'TC6','FWTC_C_Avg(7)':'TC7',}, inplace=True)
ATMe_cali.loc[:,'DOY'] = pd.Series(ATMe_cali.index.dayofyear, index=ATMe_cali.index)


#calibrate all Tc against Ts and put them in ATMs_cali dataframe
#TC4 and TC5 are broken and do not need to be 
coeffeTC=[]

for tce,tcs,ixs in zip(ATMe.columns[[1,2,3,6,7]],ATMs.columns[[1,2,3,6,7]],[0,1,2,4,5]):  #the first two idex with "names" and the last (ixs) with indeces
    coefe,rre  = poly.polyfit(ATMe[tce], ECe.Ts, 1,full=True)  #coeff, r  with r= [sumofresiduals, rank, singular_values, rcond] : list
    coeffeTC.append(np.append(tce+'end',coefe).reshape((1,3)))
    ATMe_cali[tce] = poly.polyval(ATMe[tce], coefe)
    if plotting:
        fig,ax = plt.subplots(2,1,figsize=(8,12))
        ax[0].plot(ATMe[tce],ECe.Ts,'o',color='goldenrod',markersize=5) #30cm
        ax[0].plot(ATMe[tce], ATMe_cali[tce],color='navy')
        ax[0].set_xlabel('TC temp C')
        ax[0].set_ylabel('T sonic temp C')
        ax[0].set_title(tce+'end')
        ax[0].legend(['data','fit'])
        ax[0].text(-20,-10,('coefs:  m='+str(coefe[1])+' c='+str(coefe[0])))
        
        # plot in subplot 1
        ax02= inset_axes(ax[0],height="25%", # set height
                              width="25%", # and width
                              loc=4)
        ax02.hist(ECe.Ts-ATMe_cali[tce])
        ax02.set_xlabel=('count')
        ax02.set_title('Residuals  sum= '+str(rr[0]))
        
    
        
        #subplot 2
        ax[1].plot(ATMs[tcs],ECs.Ts,'o',color='goldenrod',markersize=5) #30cm
        ax[1].plot(ATMs[tcs], ATMs_cali[tcs],color='navy')
        ax[1].set_xlabel('TC temp C')
        ax[1].set_ylabel('T sonic temp C')
        ax[1].set_title(coeffsTC[ixs][0,0])
        ax[1].legend(['data','fit'])
        ax[1].text(-32,-22,('coefs:  m='+str(coeffsTC[ixs][0,2])+' c='+str(coeffsTC[ixs][0,1])))
        
        #plot in subplot 2
        ax12=inset_axes(ax[1],height="25%", # set height
                              width="25%", # and width
                              loc=4)
        ax12.hist(ECs.Ts-ATMs_cali[tcs])
        ax12.set_title='Residuals  sum= '+str(rr[0])
        ax12.set_xlabel=('count')
        fig.set_tight_layout
        fig.show()


#%% combined calibration periods

ATMcombi=pd.concat((ATMs,ATMe))  #TC4 & 5 are broken; not usable 
ECcombi=pd.concat((ECs,ECe))


c_combi=np.zeros((2,5))

for ix,tc in enumerate(ATMcombi.columns[[1,2,3,6,7]]):  #counter,value
    c_combi[:,ix] = poly.polyfit(ATMcombi[tc],ECcombi.Ts,  1,full=False) #polyfit(x,y,degree)

#to check if calibration is visually okay
TCcalibrated_combi=pd.DataFrame(data=np.zeros_like(ATMcombi[ATMcombi.columns[[1,2,3,6,7]]]),index=ATMcombi.index,columns=ATMcombi.columns[[1,2,3,6,7]])  # 1st row as the column names

#load full season 10 min data
path_fullseason='/Volumes/U/EC_DATA/20180803_endofmeasurements/'
ATM_fullseason=loadperiod(path_fullseason,'CR3000_Table_ae.dat',"2018-05-19 11:50:00","2018-08-03 11:00:00")
ATM_fullseason.rename(columns={'RECORD': 'DOY','FWTC_C_Avg(1)':'TC1','FWTC_C_Avg(2)':'TC2','FWTC_C_Avg(3)':'TC3','FWTC_C_Avg(4)':'TC4','FWTC_C_Avg(5)':'TC5','FWTC_C_Avg(6)':'TC6','FWTC_C_Avg(7)':'TC7',}, inplace=True)

#apply calibration to full season TC record which is named ATM_fullseason and save it to dataframe TCcalibrated_all which can be exported
TCcalibrated_all=pd.DataFrame(data=ATM_fullseason[ATMcombi.columns[[0,1,2,3,6,7]]],index=ATM_fullseason.index,columns=ATMcombi.columns[[0,1,2,3,6,7]])  # 1st row as the column names


for ix,tc in enumerate(ATMcombi.columns[[1,2,3,6,7]]):  #counter,value
    TCcalibrated_combi[tc] = poly.polyval(ATMcombi[tc], c_combi[:,ix])  #eigetnliche calibrierung
    TCcalibrated_all[tc] = poly.polyval(ATM_fullseason[tc], c_combi[:,ix])  #eigetnliche calibrierung
    if plotting:
        fig,ax1 = plt.subplots(figsize=(16,8))
        ax1.plot(ATMcombi[tc],ECcombi.Ts,'o',color='goldenrod',markersize=5) 
        ax1.plot(ATMcombi[tc], TCcalibrated_combi[tc],color='navy')
        ax1.set_xlabel('TC temp C')
        ax1.set_ylabel('T sonic temp C')
        ax1.set_title(tc+' combi')
        ax1.legend(['data','fit'])
            
        ax2=fig.add_axes([0.6,0.15,0.25,0.25],title='Residuals      sum= '+str(rr[0]))
        ax2.hist(ECcombi.Ts-TCcalibrated_combi[tc])
        ax2.set_xlabel=('count')
        plt.figtext(0.2,0.65,('coefs:  m='+str(c_combi[1,ix])+' c='+str(c_combi[0,ix])))
        fig.set_tight_layout
        fig.show()

#export calibrated 10 min data
exportdf(TCcalibrated_all,'/Users/swa048/Documents/Data/EGRIP_Data/EGRIP_TC/','TC_EGRIP18_calibrated')