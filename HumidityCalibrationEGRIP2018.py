#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 31 11:51:37 2020
Humidity calibration for EGRIp 2018 season

26.05.2021 but I use a copy of the 2019 humidiy calibration script 

we calibrate the new 2140 Picarro (averaged to 1h) in the period '2018-06-10 12:00' 2018-07-01 21:00' where it was only measuring at 2m level
to the geus data from 2m.

The 2120 needs its own calibration 

2140 broken also needs its own .. :/ or if I do not use it for the snowiso model maybe I do not need it


@author: swa048
"""

import numpy as np
import pandas as pd
from scipy import signal
from scipy import stats as scistats
import datetime as dt
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from loadandconvert import loadEC,loadECperiod,removeoutliers,dfperiodfun,deltatogas,loadfiledat,plot1,exportdf
import matplotlib.cm as cm
#import windrose
import numpy.polynomial.polynomial as poly

import typhon.physics as ty

#interpolation to a new datetimeindex
def interp(df, new_index):
    """Return a new DataFrame with all columns values interpolated
    to the new_index values."""
    df_out = pd.DataFrame(index=new_index)
    df_out.index.name = df.index.name

    for colname, col in df.iteritems():
        df_out[colname] = np.interp(new_index, df.index, col)

    return df_out

# or this method
#upsample to 30 min interval and then mean daily
#new_index=pd.DatetimeIndex(start=dfPIC.index[0].strftime('%Y-%m-%d %H'),end=dfPIC.index[-1].strftime('%Y-%m-%d %H'),freq='30Min',yearfirst=True)
#dummy_frame = pd.DataFrame(np.NaN, index=new_index, columns=dfPIC.columns)
#dfPIC30=dfPIC.combine_first(dummy_frame).interpolate('time',limit=3).resample('30Min').asfreq() #this interpolated the raw data to the dummy_frame index values, and then resamples to the dummy_frame values
# df.combine_first(dummy_frame) adds the evenly spaced timeindex (full seconds) to the dataframe, 
# then interpolate, interpolates it to all of those times (dataframe is now doubled in length) but leave out gaps that are larger than 30s (60 cells) long 
#and resample cuts it down to 1s- resolution again.



#%% new humidity calibration for Vapour 2018


''' !!! small note:
    HC thinks he is calibrating against mixing ratio but he is in fact calibrating against specific humidity!. minor difference
    Picarro is calibrated against the geus humidity within the period where Picarro was only measuring at 2m level
    Picarro is calibrated against the specific humidity measuremens as it is also converted to specific humidity  using HCs formula
    in the end, it does not really matter if its specific huidity or mixing ratio'''

# import geus
# to get HC's calibration values, you need to import the old GEUS datafile /Users/swa048/forServer/Meteo/PROMICE_GEUS/2019/EGP_2019_v03/EGP_hour_v03.txt
# to get the values I have been using, you need the updated GEUS datagile
mydateparser = lambda x: dt.datetime.strptime(x, "%Y %m %d %H")
geus_all=pd.read_csv('/Users/swa048/forServer/Meteo/PROMICE_GEUS/2019/EGP_2020_06_24_v03/EGP_hour_v03.txt',delim_whitespace=True,index_col=(0),parse_dates=[[0, 1,2,3]],infer_datetime_format=True ,date_parser=mydateparser,na_values=['NAN',-999.00],usecols=(0,1,2,3,4,6,7,8,9,10,11,12,13,14,24,25))
geus_all.rename(inplace=True,columns={ 
        'AirPressure(hPa)': 'p',
        'AirTemperature(C)': 'temp',
        'RelativeHumidity(%)':'relhum', #is calculated with respect to ice
        'WindSpeed(m/s)':'wsp',
        'SpecificHumidity(g/kg)':'spechum', #g/kg ambient air calculated from RH (not dry air but moist air)
        'HeightSensorBoom(m)':'ranger_high',
        'HeightStakes(m)':'ranger_low',
        'LatentHeatFlux(W/m2)':'LE',
        'SensibleHeatFlux(W/m2)':'H',
        'WindDirection(d)':'wd'})

geus_h=geus_all.loc[pd.to_datetime('2018-06-10 12:00'):pd.to_datetime('2018-07-01 21:00'),:].copy() #thats the period when Pic was only running at 2m level


B=621.9907 #g/kg
''''# This is HC's convertion approach, this is actually used in all paper 1 calculations. It's almost the same as the laternative calculation below
geus_h.loc[:,'Pw']= (geus_h.spechum/B)*(geus_h.p*100)/(1+(geus_h.spechum/B)) #Water Vapour Pressure Pa from Vaisala eq 14, I assume spechum is the same as mixing ratio

C=2.16679 #gK/J
geus_h.loc[:,'abshum']= C*geus_h.Pw/(geus_h.temp+273.15) #absolute humidity in g/m3 from Vaisala eq 17, assuming ideal gas law
'''
# alternative way of calculating as of 15.02.2021 (is more or less the same as above, except that you do not use eq14 with assuming spec humidity is mixing ratio)
geus_h.loc[:,'abshum']= geus_h.spechum*ty.density(geus_h.p*100, T=geus_h.temp+273.15)  #calculating density assuming ideal gas law and dry air constant 


# see if what I do is comparable to the typhon absolute huidity calculation --> they dont have that conversion, even better :D

# see what the difference between specific humidity and mixing ratio is in g/kg
geus_h.loc[:,'mr']=ty.specific_humidity2mixing_ratio(geus_h.spechum/1000)*1000


#%% PIC 2140

# PIC, it seems HC is actially calculating the specific humidity with his calculation (instead of mixing ratio as he thinks), HOWEVER, it does not make a difference
dfPIC2m10min=pd.read_csv('/Volumes/gfi_snowiso/Vapour/2018/processed_data_combined/EG18_level180cm.txt',index_col=0,parse_dates=True,na_values=['NAN'])
dfPIC2m10min=dfPIC2m10min.loc[pd.to_datetime('2018-06-11 12:00'):pd.to_datetime('2018-07-01 21:00'),:] #only use specific period where the new 2140 was measuring at 2m only
dfPIC2m10min.loc[:,'shum']=dfPIC2m10min.H2O*1e-6*B #you could also use deltatogas function, dfPIC2m.mr is g/kg air --> specific humidity


#interpolate to 1h mean
dfPIC2m=dfPIC2m10min.resample('1h').mean()



# CHECK what typhon does --> Conclusion: What I have done is pretty close to what typhon would do, so no need to change my code again
# see what typhon makes from the volume mixing ratio (ppmv) 
dfPIC2m.loc[:,'ty_mr']=ty.vmr2mixing_ratio(dfPIC2m.H2O*1e-6)*1e3 # to get g/kg

dfPIC2m.loc[:,'ty_mr_from_mr']=ty.specific_humidity2mixing_ratio(dfPIC2m.shum/1000)*1000  #assume HC formula is actually for spec humidity and not for mr
dfPIC2m.loc[:,'ty_spechum']=ty.vmr2specific_humidity(dfPIC2m.H2O*1e-6)*1e3

# combine to get rid of nan
dfcombiPIC=pd.concat([dfPIC2m.shum,geus_h.spechum],axis=1).dropna() # dfcombiPIC dfPIC2m.shum in g/kg and dfgeus.spechum in g/kg 


#%% calibration of Picarro 2140 EASTGRIP 2019

plt.figure()
plt.plot(dfPIC2m10min.shum,label='PIC original',Marker='s',LineStyle='-',color='darkred',MarkerSize=2)
plt.plot(dfPIC2m.shum,label='PIC resampled',Marker='*',LineStyle=':',color='blue',MarkerSize=2)
plt.plot(geus_h.spechum,label='geus',Marker='s',LineStyle='--',color='orange',MarkerSize=2)
plt.legend()

plt.figure()
plt.scatter(dfcombiPIC.shum,dfcombiPIC.spechum,s=12,color='darkred') #shum = PIC, spechum=geus
plt.xlabel('PIC g/kg air')
plt.ylabel('geus g/kg air')
plt.title('specific humidity calibration PIC 2140')



coefPIC,rPIC = poly.polyfit(dfcombiPIC.shum,dfcombiPIC.spechum,deg=1,full=True) #polyfit(x,y,degree) .shum = PIC, spechum=geus
dfPIC2m.loc[:,'qcali'] = poly.polyval(dfPIC2m.shum, coefPIC)  #eigetnliche calibrierung


fig2,ax=plt.subplots(1,1)
ax.scatter(dfcombiPIC.shum,dfcombiPIC.spechum,color='darkred',s=15,label='original')
ax.plot(dfPIC2m.shum,dfPIC2m.qcali,label='calibration',color='greenyellow')
ax.text(2,1,('coefs:  m='+str(np.round(coefPIC[1],4))+' c='+str(np.round(coefPIC[0],4))))
ax.set_ylabel('geus g/kg (from rel.hum)')
ax.set_xlabel('PIC 2m g/kg')
ax.set_title('specific humidity PIC 2140')
fig2.savefig('/Users/swa048/forServer/Vapour/2018/humcal_PIC2140geusEGRIP2018.png')



fig,ax1,ax2=plot1()
ax1.plot(dfPIC2m10min.shum,label='PIC original')
ax1.plot(dfPIC2m.qcali,label='PIC calibrated')
ax1.plot(geus_h.spechum,label='geus data')
ax1.legend()
fig.autofmt_xdate()
ax1.set_ylabel('specific humidity g/kg ')
fig.savefig('/Users/swa048/forServer/Vapour/2018/timeseries_PIC2140geusEGRIP2018.png')


# lin fit forced through zero --> we need the normal fit 
# Our model is y = a * x, so things are quite simple, in this case...
# x needs to be a column vector instead of a 1D vector for this, however.
#x=dfcombiIrga.q
#y=dfcombiIrga.abshum
#x = x[:,np.newaxis]
#a, _, _, _ = np.linalg.lstsq(x, y)

#plt.figure()
#plt.plot(x, y, 'bo')
#plt.plot(x, a*x, 'r-')
#plt.show()
#%% PIC 2120

#here we interpolate the geus to the times of Picarro because I can not resample only 1 datapoint to 1 hour

# PIC 2120 #PIC 2120 data from 24.05.2018 17:16 - 10.06.2018 12:51
dfPIC2m10min2=pd.read_csv('/Volumes/gfi_snowiso/Vapour/2018/processed_data_combined/EG18_level180cm.txt',index_col=0,parse_dates=True,na_values=['NAN'])
dfPIC2120=dfPIC2m10min2.loc[pd.to_datetime('2018-05-24 17:16'):pd.to_datetime('2018-06-10 13:00'),:].copy() #only use specific period where the 2120 was measuring at 2m only
dfPIC2120.loc[:,'shum']=dfPIC2120.H2O*1e-6*B #you could also use deltatogas function, dfPIC2m.mr is g/kg air --> specific humidity

geus_h=geus_all.loc[pd.to_datetime('2018-05-24 17:16'):pd.to_datetime('2018-06-10 13:00'),:].copy() #thats the period when Pic was only running at 2m level

#%% interpoalte

# geus to new timevec with Picarro spacing 
geus_new=interp(geus_h,dfPIC2120.index)

# combine to get rid of nan
dfcombiPIC2120=pd.concat([dfPIC2120.shum,geus_new.spechum],axis=1).dropna() # dfcombiPIC dfPIC2m.shum in g/kg and dfgeus.spechum in g/kg 

#%% calibration of Picarro 2120 EASTGRIP 2018

plt.figure()
plt.plot(dfPIC2120.shum,label='PIC original',Marker='s',LineStyle='-',color='darkred',MarkerSize=2)
plt.plot(geus_h.spechum,label='geus',Marker='s',LineStyle='--',color='orange',MarkerSize=2)
plt.plot(geus_new.spechum,label='geus resampled',Marker='*',LineStyle=':',color='blue',MarkerSize=2)
plt.legend()

plt.figure()
plt.scatter(dfcombiPIC2120.shum,dfcombiPIC2120.spechum,s=12,color='darkred') #shum = PIC, spechum=geus
plt.xlabel('PIC g/kg air')
plt.ylabel('geus g/kg air')
plt.title('specific humidity calibration PIC 2120')



coefPIC2120,rPIC2120 = poly.polyfit(dfcombiPIC2120.shum,dfcombiPIC2120.spechum,deg=1,full=True) #polyfit(x,y,degree) .shum = PIC, spechum=geus
dfPIC2120.loc[:,'qcali'] = poly.polyval(dfPIC2120.shum, coefPIC2120)  #eigetnliche calibrierung


fig2,ax=plt.subplots(1,1)
ax.scatter(dfcombiPIC2120.shum,dfcombiPIC2120.spechum,color='darkred',s=15,label='original')
ax.plot(dfPIC2120.shum,dfPIC2120.qcali,label='calibration',color='greenyellow')
ax.text(2,1,('coefs:  m='+str(np.round(coefPIC2120[1],4))+' c='+str(np.round(coefPIC2120[0],4))))
ax.set_ylabel('geus g/kg (from rel.hum)')
ax.set_xlabel('PIC 2m g/kg')
ax.set_title('specific humidity PIC 2120')
fig2.savefig('/Users/swa048/forServer/Vapour/2018/humcal_PIC2120geusEGRIP2018.png')



fig,ax1,ax2=plot1()
ax1.plot(dfPIC2120.shum,label='PIC original')
ax1.plot(dfPIC2120.qcali,label='PIC calibrated')
ax1.plot(geus_h.spechum,label='geus data')
ax1.set_ylabel('specific humidity g/kg ')
ax1.legend()
fig.autofmt_xdate()
fig.savefig('/Users/swa048/forServer/Vapour/2018/timeseries_PICb2120EGRIP2018.png')
#%% PIC 2140 broken

#here we interpolate the geus to the times of Picarro because I can not resample only 1 datapoint to 1 hour

# PIC 2140 ## broken PIC data from 11.05. 00:00 - 23.05. 18:08
dfPICbroken=dfPIC2m10min2.loc[pd.to_datetime('2018-05-11 00:00'):pd.to_datetime('2018-05-23 18:08'),:].copy() #only use specific period where the 2120 was measuring at 2m only
dfPICbroken.loc[:,'shum']=dfPICbroken.H2O*1e-6*B #you could also use deltatogas function, dfPIC2m.mr is g/kg air --> specific humidity

geus_h2=geus_all.loc[pd.to_datetime('2018-05-11 00:00'):pd.to_datetime('2018-05-23 18:08'),:].copy() #thats the period when Pic was only running at 2m level

# geus to new timevec with Picarro spacing 
geus_new=interp(geus_h2,dfPICbroken.index)

# combine to get rid of nan
dfcombiPICbroken=pd.concat([dfPICbroken.shum,geus_new.spechum],axis=1).dropna() # dfcombiPIC dfPIC2m.shum in g/kg and dfgeus.spechum in g/kg 

#%% calibration of Picarro 2140 broken EASTGRIP 2018

plt.figure()
plt.plot(dfPICbroken.shum,label='PIC original',Marker='s',LineStyle='-',color='darkred',MarkerSize=2)
plt.plot(geus_h2.spechum,label='geus',Marker='s',LineStyle='--',color='orange',MarkerSize=2)
plt.plot(geus_new.spechum,label='geus resampled',Marker='*',LineStyle=':',color='blue',MarkerSize=2)
plt.legend()

plt.figure()
plt.scatter(dfcombiPICbroken.shum,dfcombiPICbroken.spechum,s=12,color='darkred') #shum = PIC, spechum=geus
plt.xlabel('PIC g/kg air')
plt.ylabel('geus g/kg air')
plt.title('specific humidity calibration PIC 2140 broken')



coefPICbroken,rPICbroken = poly.polyfit(dfcombiPICbroken.shum,dfcombiPICbroken.spechum,deg=1,full=True) #polyfit(x,y,degree) .shum = PIC, spechum=geus
dfPICbroken.loc[:,'qcali'] = poly.polyval(dfPICbroken.shum, coefPICbroken)  #eigetnliche calibrierung


fig2,ax=plt.subplots(1,1)
ax.scatter(dfcombiPICbroken.shum,dfcombiPICbroken.spechum,color='darkred',s=15,label='original')
ax.plot(dfPICbroken.shum,dfPICbroken.qcali,label='calibration',color='greenyellow')
ax.text(1.25,0.4,('coefs:  m='+str(np.round(coefPICbroken[1],4))+' c='+str(np.round(coefPICbroken[0],4))))
ax.set_ylabel('geus g/kg (from rel.hum)')
ax.set_xlabel('PIC 2m g/kg')
ax.set_title('specific humidity PIC 2140 broken')
fig2.savefig('/Users/swa048/forServer/Vapour/2018/humcal_PICbrokengeusEGRIP2018.png')



fig,ax1,ax2=plot1()
ax1.plot(dfPICbroken.shum,label='PIC original')
ax1.plot(dfPICbroken.qcali,label='PIC calibrated')
ax1.plot(geus_h2.spechum,label='geus data')
ax1.legend()
fig.autofmt_xdate()
ax1.set_ylabel('specific humidity g/kg ')
fig.savefig('/Users/swa048/forServer/Vapour/2018/timeseries_PICbrokengeusEGRIP2018.png')

#%% export humidity calibrated dataframes

# 1.8m data vapor
df2mraw=pd.read_csv('/Users/swa048/forServer/Vapour/2018/processed_data_combined/EG18_level180cm.txt',index_col=0,parse_dates=True,na_values=['NAN'])
df2mraw.rename(inplace=True,columns={ 'H2O' : 'q'})
df2m2120=deltatogas(df2mraw.loc[pd.to_datetime('2018-05-24 17:16'):pd.to_datetime('2018-06-10 13:00'),:],m=0.9294, c=-0.062)
df2m2140=deltatogas(df2mraw.loc[pd.to_datetime('2018-06-11 17:45'):pd.to_datetime('2018-08-03 14:23'),:],m=0.7871, c=-0.0631)
df2mbroken=deltatogas(df2mraw.loc[pd.to_datetime('2018-05-11 00:00'):pd.to_datetime('2018-05-23 18:08'),:],m=0.7804, c=-0.0021)
df2m=pd.concat([df2mbroken,df2m2120,df2m2140],axis=0)

# 7m data vapor
df7mraw=pd.read_csv('/Users/swa048/forServer/Vapour/2018/processed_data_combined/EG18_level7m.txt',index_col=0,parse_dates=True,na_values=['NAN'])
df7mraw.rename(inplace=True,columns={ 'H2O' : 'q'})
df7m2120=deltatogas(df7mraw.loc[pd.to_datetime('2018-05-24 17:16'):pd.to_datetime('2018-06-10 13:00'),:],m=0.9294, c=-0.062)
df7m2140=deltatogas(df7mraw.loc[pd.to_datetime('2018-06-11 17:45'):pd.to_datetime('2018-08-03 14:23'),:],m=0.7871, c=-0.0631)
df7mbroken=deltatogas(df7mraw.loc[pd.to_datetime('2018-05-11 00:00'):pd.to_datetime('2018-05-23 18:08'),:],m=0.7804, c=-0.0021)
df7m=pd.concat([df7mbroken,df7m2120,df7m2140],axis=0)

# 80cm data vapor
df80cmraw=pd.read_csv('/Users/swa048/forServer/Vapour/2018/processed_data_combined/EG18_level80cm.txt',index_col=0,parse_dates=True,na_values=['NAN'])
df80cmraw.rename(inplace=True,columns={ 'H2O' : 'q'})
df80cm2120=deltatogas(df80cmraw.loc[pd.to_datetime('2018-05-24 17:16'):pd.to_datetime('2018-06-10 13:00'),:],m=0.9294, c=-0.062)
df80cm2140=deltatogas(df80cmraw.loc[pd.to_datetime('2018-06-11 17:45'):pd.to_datetime('2018-08-03 14:23'),:],m=0.7871, c=-0.0631)
df80cmbroken=deltatogas(df80cmraw.loc[pd.to_datetime('2018-05-11 00:00'):pd.to_datetime('2018-05-23 18:08'),:],m=0.7804, c=-0.0021)
df80cm=pd.concat([df80cmbroken,df80cm2120,df80cm2140],axis=0)

# 30cm data vapor
df30cmraw=pd.read_csv('/Users/swa048/forServer/Vapour/2018/processed_data_combined/EG18_level30cm.txt',index_col=0,parse_dates=True,na_values=['NAN'])
df30cmraw.rename(inplace=True,columns={ 'H2O' : 'q'})
df30cm2120=deltatogas(df30cmraw.loc[pd.to_datetime('2018-05-24 17:16'):pd.to_datetime('2018-06-10 13:00'),:],m=0.9294, c=-0.062)
df30cm2140=deltatogas(df30cmraw.loc[pd.to_datetime('2018-06-11 17:45'):pd.to_datetime('2018-08-03 14:23'),:],m=0.7871, c=-0.0631)
df30cmbroken=deltatogas(df30cmraw.loc[pd.to_datetime('2018-05-11 00:00'):pd.to_datetime('2018-05-23 18:08'),:],m=0.7804, c=-0.0021)
df30cm=pd.concat([df30cmbroken,df30cm2120,df30cm2140],axis=0)
#%% look at the calibrated 


fig,ax1,ax2=plot1()
ax1.plot(geus_all.loc[pd.to_datetime('2018-05-11 00:00'):pd.to_datetime('2018-08-03 14:23'),'spechum'],label='geus',c='k') 
ax1.plot(df2m.loc[pd.to_datetime('2018-05-11 00:00'):pd.to_datetime('2018-05-23 18:08'),'q'],label='calibrated 2140 broken')
ax1.plot(df2m.loc[pd.to_datetime('2018-05-24 17:16'):pd.to_datetime('2018-06-10 13:00'),'q'],label='calibrated 2120')
ax1.plot(df2m.loc[pd.to_datetime('2018-06-11 17:45'):pd.to_datetime('2018-08-03 14:23'),'q'],label='calibrated 2140')
ax1.legend()
ax1.set_title('180cm level ')
ax1.set_ylabel('specific humidity [g/kg]',fontsize=18)
fig.autofmt_xdate()
fig.savefig('/Users/swa048/forServer/Vapour/2018/timeseries_allPIChumcalibrated_2m.png')

#%% #those are the dataframes that have humidity as specific humidity and are calibrated.. they rely on the calibration further down in this script
exportdf(df7m,'/Users/swa048/forServer/Vapour/2018/processed_data_combined/','EG18_humcal7m')
exportdf(df30cm,'/Users/swa048/forServer/Vapour/2018/processed_data_combined/','EG18_humcal30cm')
exportdf(df80cm,'/Users/swa048/forServer/Vapour/2018/processed_data_combined/','EG18_humcal80cm')
exportdf(df2m,'/Users/swa048/forServer/Vapour/2018/processed_data_combined/','EG18_humcal180cm')
