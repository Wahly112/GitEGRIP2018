# -*- coding: utf-8 -*-
"""
Created on Mon Apr 30 13:02:33 2018

@author: Sonja

plot daily cycles of 
temperature
fluxes
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import pandas as pd
import numpy.polynomial.polynomial as poly
import datetime as dt

# =============================================================================
# plt.style.use('ggplot')
# plt.rc('pgf',  texsystem='pdflatex')
# plt.rc('text', usetex=True)
# =============================================================================
plt.rcParams['text.latex.preamble']=[r"\usepackage{wasysym}"]



#from calibration_EGRIP2018 import coeffsTC

coeffsTC=[np.array([['TC1start', '-1.57362075501', '0.959301058023']],
      dtype='|S32'), np.array([['TC2start', '-1.40215716179', '0.965248243814']],
      dtype='|S32'), np.array([['TC3start', '-1.30694164465', '0.968106549167']],
      dtype='|S32'), np.array([['TC5start', '-1.01507974412', '0.97468185986']],
      dtype='|S32'), np.array([['TC6start', '-1.1588187557', '0.965786277967']],
      dtype='|S32'), np.array([['TC7start', '-1.54208009458', '0.955341455148']],
      dtype='|S32'), np.array([['TC1end', '-1.57362075501', '0.959301058023']],
      dtype='|S32'), np.array([['TC2end', '-1.40215716179', '0.965248243814']],
      dtype='|S32'), np.array([['TC3end', '-1.30694164465', '0.968106549167']],
      dtype='|S32'), np.array([['TC6end', '-1.1588187557', '0.965786277967']],
      dtype='|S32'), np.array([['TC7end', '-1.54208009458', '0.955341455148']],
      dtype='|S32')]


#%%
path='/Volumes/U/EC_DATA/20180803_endofmeasurements/'
path_plot='/Volumes/U/EC_processed/plots/daily/'
filename='CR3000_ec_scfd.dat'
ATMfilename='CR3000_Table_ae.dat'

np.warnings.filterwarnings('ignore')  # when comparing to nan tha warning "invalid value encountered" will aprea. This way it is only given once

plt.close('all')
plotting = 0

''' I think this is unnecessary and from the early beginning of my Python coding'''
##read timestamp from datafile and remove quotation marks ""
#date_string=np.genfromtxt(path+filename,dtype=None, delimiter=',', skip_header=4, skip_footer=0, usecols =(0))
#for i in range(len(date_string)):
#     date_string[i]=date_string[i].replace('"', '')  #remove the quotation marks
#date=np.array(date_string, dtype='datetime64')
##use (a-b)/np.timedelta64(1,'h') to get difference in h,s,m...
#
## read variable names from file and save it in var, then load data and access through names
#file = open(path+filename, 'r') #read only
#var = file.readlines()[1] #names for columns
#file.close()
#var=var.replace('"','')
#var=var.replace('\n','')
#var= var.split(',')
#data=np.genfromtxt(path+filename,delimiter=',',names=var, skip_header=4, skip_footer=0)
##access array by data['LE'] or data[var[3]]

#and import 10 min averages as dataframe
df=pd.read_csv(path+filename,index_col=0,skiprows=[0,2,3],parse_dates=True,na_values=['NAN',-7999],dtype={'RECORD': np.int64, 'LE': np.float64, 'Hs': np.float64, 'tau': np.float64, 'u_star': np.float64, 'Ux_Ux': np.float64, 'Ux_Uy': np.float64, 'Ux_Uz': np.float64, 'Ux_ln_vh': np.float64, 'Ux_Ts': np.float64, 'Uy_Uy': np.float64, 'Uy_Uz': np.float64, 'Uy_ln_vh': np.float64, 'Uy_Ts': np.float64, 'Uz_Uz': np.float64, 'Uz_ln_vh': np.float64, 'Uz_Ts': np.float64, 'ln_vh_ln_vh': np.float64, 'ln_vh_Ts': np.float64, 'Ts_Ts': np.float64, 'Ux_Avg': np.float64, 'Uy_Avg': np.float64, 'Uz_Avg': np.float64, 'ln_vh_Avg': np.float64, 'Ts_Avg': np.float64, 'horiz_wind_spd': np.float64, 'result_wind_spd': np.float64, 'wind_dir': np.float64, 'wind_dir_std_dev': np.float64, 'wnd_dir_compass': np.float64, 'vh_Avg': np.float64, 'n_Tot': np.float64, 'del_T_f_Tot': np.float64, 'track_f_Tot': np.float64, 'amp_h_f_Tot': np.float64, 'amp_l_f_Tot': np.float64} )
df.rename(columns={'RECORD': 'DOY'}, inplace=True)
t0=np.datetime64(str(df.index[0])[0:4]+'-01-01')
df.loc[:,'DOY'] = pd.Series((df.index-t0)/dt.timedelta(days=1), index=df.index) 
#df.rename(columns={'ln_vh_Avg': 'hum'}, inplace=True)
#df.loc[:,'hum']=pd.Series(df.hum/(-0.105),index=df.index)





#import atmospheric parameters 
atm=pd.read_csv(path+ATMfilename,index_col=0,skiprows=[0,2,3],parse_dates=True,na_values=['NAN',-7999],dtype={'RECORD': np.int64, 'FWTC_C_Avg(1)': np.float64, 'FWTC_C_Avg(2)': np.float64, 'FWTC_C_Avg(3)': np.float64, 'FWTC_C_Avg(4)': np.float64, 'FWTC_C_Avg(5)': np.float64, 'FWTC_C_Avg(6)': np.float64, 'FWTC_C_Avg(7)': np.float64, 'WS_ms_Avg(1)': np.float64, 'WS_ms_Avg(2)': np.float64, 'WS_ms_Avg(3)': np.float64, 'WS_ms_Max(1)': np.float64, 'WS_ms_Max(2)': np.float64, 'WS_ms_Max(3)': np.float64, 'BattV': np.float64, 'PTemp_C': np.float64} )
atm.rename(columns={'RECORD': 'DOY','FWTC_C_Avg(1)':'TC1','FWTC_C_Avg(2)':'TC2','FWTC_C_Avg(3)':'TC3','FWTC_C_Avg(4)':'TC4','FWTC_C_Avg(5)':'TC5','FWTC_C_Avg(6)':'TC6','FWTC_C_Avg(7)':'TC7',}, inplace=True)
atm.loc[:,'DOY'] = pd.Series(atm.index.dayofyear, index=atm.index)


#calibrate with values from end where possible
for tc in atm.columns[[1,2,3,5,6,7]]:
    if tc=='TC1':
        c=coeffsTC[6][0,1:3].astype(np.float)
        atm[tc] = poly.polyval(atm[tc], c)
    if tc=='TC2':
        c=coeffsTC[7][0,1:3].astype(np.float)
        atm[tc] = poly.polyval(atm[tc], c)
    if tc=='TC3':
        c=coeffsTC[8][0,1:3].astype(np.float)
        atm[tc] = poly.polyval(atm[tc], c)
    if tc=='TC5':
        c=coeffsTC[3][0,1:3].astype(np.float)
        atm[tc] = poly.polyval(atm[tc], c)
    if tc=='TC6':
        c=coeffsTC[9][0,1:3].astype(np.float)
        atm[tc] = poly.polyval(atm[tc], c)    
    if tc=='TC7':
        c=coeffsTC[10][0,1:3].astype(np.float)
        atm[tc] = poly.polyval(atm[tc], c)     



#%% contour plot 

dateRange = pd.date_range('05/20/2018', '08/03/2018', freq='D')
dateRange=dateRange.format(formatter=lambda x: x.strftime('%Y-%m-%d'))



for day in dateRange:
    daycp=atm[day]
    daycp.index=daycp.index.time
    daycp.index=daycp.index.format(formatter=lambda x: x.strftime('%H:%M'))
    xlist=daycp.index  
    ylist=np.array((2,5,10,30,60,180))
    X,Y=np.meshgrid(xlist,ylist)
    ZT=daycp[daycp.columns[[5,6,7,1,2,3]]]
    Z=np.transpose(ZT)
    fig=plt.figure(figsize=(16,10))
    cp=plt.contourf(X,Y,Z,14,cmap='RdBu_r')
    locs,labels=plt.xticks()
    plt.xticks(np.arange(0,144,6),pd.date_range(daycp.index[0],daycp.index[-1],freq='H').format(formatter=lambda x: x.strftime('%H:%M')),rotation='vertical')
    cbar=plt.colorbar(cp)
    cbar.ax.set_title(u'T in $^\circ$C')
    plt.xlabel('UTC')
    plt.ylabel('heightabove ground (cm)')
    plt.title(day+' / DOY: '+str(daycp.DOY[0]))
    #fig.savefig(path_plot+'Tcontour_'+day+".pdf")#, bbox_inches='tight')
    plt.close()

 

if plotting:       
        
    dateRange = pd.date_range('05/19/2018', '08/02/2018', freq='D')
    dateRange=dateRange.format(formatter=lambda x: x.strftime('%Y-%m-%d'))
    for day in dateRange:
        dayzzz=atm[day]
        fig=plt.figure(figsize=(16,10))
        dayzzz.TC5.plot(color='blueviolet') #2cm
        dayzzz.TC6.plot(color='crimson') #5cm
        dayzzz.TC7.plot(color='indianred') #10cm
        dayzzz.TC1.plot(color='goldenrod') #30cm
        dayzzz.TC2.plot(color='yellowgreen') #60cm
        dayzzz.TC3.plot(color='darkcyan') #180cm
        #dayzzz.TC4.plot(color='darkcyan')
        plt.title(day+' / DOY: '+str(dayzzz.DOY[0]))
        plt.legend(['2cm','5cm','10cm','30cm','60cm','180cm'])
        plt.ylabel('degrees C')
        plt.xlabel('')
        #fig.savefig(path_plot+'Tgrad_'+day+".pdf")#, bbox_inches='tight')
        plt.close()









    dateRange = pd.date_range('05/19/2018', '08/03/2018', freq='D')
    dateRange=dateRange.format(formatter=lambda x: x.strftime('%Y-%m-%d'))
    for day in dateRange:
        dayxxx=df[day]
        fig=plt.figure(figsize=(8,12))
        plt.subplot(4,1,1)
        dayxxx.Ts_Avg.plot(color='crimson')
        plt.title(day+' / DOY: '+str(dayxxx.DOY[0]))
        plt.legend(['sonic Temperature'])
        plt.ylabel('degrees C')
        plt.xlabel('')
        plt.subplot(4,1,2)
        dayxxx.ln_vh_Avg.plot(ax=plt.gca(),color='royalblue')
        plt.legend(['KH20 humidity'])
        plt.ylabel('ln vh')
        plt.xlabel('')
        plt.subplot(4,1,3)
        dayxxx.LE.plot(ax=plt.gca(),color='olivedrab')
        plt.legend(['LE'])
        plt.ylabel('kg/(m2*s)')
        plt.xlabel('')
        plt.subplot(4,1,4)
        dayxxx.Hs.plot(ax=plt.gca(),color='coral')
        plt.legend(['Hs'])
        plt.ylabel('J/(m2s)')
        plt.tight_layout()
        #fig.savefig(path_plot+'Summary_'+day+".pdf")#, bbox_inches='tight')
        plt.close()

#OR
#load all data without the names and access as array indexing
#data=np.genfromtxt(r'C:\Users\Sonja\Documents\VAPOR4\CR3000_nottobeoverwritten\20170802\CR3000_ec_scfd.dat',delimiter=',', skip_header=4, skip_footer=0)



    

# trim LE dataset 
    df['LE'][df['LE']>50]=np.nan
    df['LE'][df['LE']<-50]=np.nan

    plt.figure() 
    plt.plot(df.index,df['LE'])
    plt.xlabel('time')
    plt.ylabel(' W/m2')
    plt.title('latent heat')
    plt.show()

    plt.figure() 
    plt.plot(df.index,df['vh_Avg'],color='steelblue')
    plt.xlabel('time')
    plt.ylabel('[mV]')
    plt.title('Kh20 10 min average')
    plt.show()


    f3=plt.figure(figsize=(18,6)) 
    plt.plot_date(df.index,df['Ts_Avg'],marker='.',color='firebrick')
    plt.gca().xaxis.set_major_locator(mdates.WeekdayLocator(byweekday=1,interval=2))
    plt.xlabel('time')
    plt.ylabel('C')
    plt.title('sonic temperature 10 min average')
    f3.autofmt_xdate()
    plt.show()


