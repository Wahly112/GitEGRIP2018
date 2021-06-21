#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 15 12:56:36 2019


# new on 26.05.2021 because wrong Date 


@author: swa048
"""


import math
import pandas as pd
import numpy as np
import datetime as dt
from scipy import interpolate
from scipy import signal
from scipy import stats as scistats
import scipy.io as sio
import scipy.io.matlab.mio5 as sio5
import matplotlib.pyplot as plt

from loadandconvert import exportdf


pathm1='/Users/swa048/forServer/Vapour/2018/PIC_broken_EGRIP/Processed_data/processed_data.mat'
pathm2='/Users/swa048/forServer/Vapour/2018/PICARRO_2120/Processed_data/processed_data_new.mat'
pathm3='/Users/swa048/forServer/Vapour/2018/Picarro_2140/Processed_data/processed_data.mat'

#import calibrated data in matrices
tdata1 = sio.loadmat(pathm1,squeeze_me=True) # broken PIC data from 11.05. 00:00 - 23.05. 18:08
tdata2= sio.loadmat(pathm2,squeeze_me=True) #PIC 2120 data from 24.05.2018 17:16 - 10.06.2018 12:51
tdata3=sio.loadmat(pathm3,squeeze_me=True) #PIC 2140 data from 11.06.2018 17:45  -  03.08. 14:23
#%% export the full dataset, not only the means


cols=['DOY','Valve','','H2O','d17O','d18O','dD']
PIC2140_2018=pd.DataFrame(data=tdata3['D_VSMOW'],columns=cols)
PIC2140_2018.index=pd.to_datetime('2018-1-1') + pd.to_timedelta(PIC2140_2018.DOY, unit='D') 
#exportdf(PIC2140_2018,'/Users/swa048/forServer/Vapour/2018/processed_data/','EGRIP2018_2140calibrated')
#%%

cols2120=['DOY','Valve','H2O','d18O','dD']
PIC2120_2018=pd.DataFrame(data=tdata2['D_VSMOW'],columns=cols2120)
PIC2120_2018.index=pd.to_datetime('2018-1-1') + pd.to_timedelta(PIC2120_2018.DOY, unit='D') 
#exportdf(PIC2120_2018,'/Users/swa048/forServer/Vapour/2018/processed_data/','EGRIP2018_2120calibrated')

#


#%%
cols2140bad=['DOY','Valve','H2O','d18O','dD']
PICbroken_2018=pd.DataFrame(data=tdata1['D_VSMOW'],columns=cols2140bad)
PICbroken_2018.index=pd.to_datetime('2018-1-1') + pd.to_timedelta(PICbroken_2018.DOY, unit='D') 
#exportdf(PICbroken_2018,'/Users/swa048/forServer/Vapour/2018/processed_data/','EGRIP2018_BROKEcalibrated')
 

#%% create dataframes according to level
cols_f=['DOY','Valve','H2O','d18O','dD','d17O','timestampPY']
lev7m=pd.DataFrame(columns=cols_f)
lev30cm=pd.DataFrame(columns=cols_f)
lev80cm=pd.DataFrame(columns=cols_f)
lev180cm=pd.DataFrame(columns=cols_f)


for mat in [tdata1,tdata2,tdata3]:
#level1:     Valve 0 = 7m
    if len(mat['mean_level1'][1,:]) <7:
        dummydf=pd.DataFrame(data=mat['mean_level1'] ,columns=['DOY','Valve','H2O','d18O','dD'])        
        lev7m=lev7m.append(dummydf,sort=True,ignore_index=True)        
    else:    
        dummydf=pd.DataFrame(data=mat['mean_level1'] ,columns=['DOY','Valve','timestampPY','H2O','d17O','d18O','dD'])        
        lev7m=lev7m.append(dummydf,sort=True,ignore_index=True)
   
#level2:     Valve 2 = 30cm
    if len(mat['mean_level2'][1,:]) <7:
        dummydf=pd.DataFrame(data=mat['mean_level2'] ,columns=['DOY','Valve','H2O','d18O','dD'])        
        lev30cm=lev30cm.append(dummydf,sort=True,ignore_index=True)        
    else:    
        dummydf=pd.DataFrame(data=mat['mean_level2'] ,columns=['DOY','Valve','timestampPY','H2O','d17O','d18O','dD'])        
        lev30cm=lev30cm.append(dummydf,sort=True,ignore_index=True)
 #level3:    Valve 4 = 80cm
    if len(mat['mean_level3'][1,:]) <7:
        dummydf=pd.DataFrame(data=mat['mean_level3'] ,columns=['DOY','Valve','H2O','d18O','dD'])        
        lev80cm=lev80cm.append(dummydf,sort=True,ignore_index=True)        
    else:    
        dummydf=pd.DataFrame(data=mat['mean_level3'] ,columns=['DOY','Valve','timestampPY','H2O','d17O','d18O','dD'])        
        lev80cm=lev80cm.append(dummydf,sort=True,ignore_index=True)
        
#level4:    Valve 8 = 180cm
    if len(mat['mean_level4'][1,:]) <7:
        dummydf=pd.DataFrame(data=mat['mean_level4'] ,columns=['DOY','Valve','H2O','d18O','dD'])        
        lev180cm=lev180cm.append(dummydf,sort=True,ignore_index=True)        
    else:    
        dummydf=pd.DataFrame(data=mat['mean_level4'] ,columns=['DOY','Valve','timestampPY','H2O','d17O','d18O','dD'])        
        lev180cm=lev180cm.append(dummydf,sort=True,ignore_index=True)
      
''' reset index with datetime '''
lev7m.timestampPY = pd.to_datetime('2018-1-1') + pd.to_timedelta(lev7m.DOY, unit='D') 
lev7m.set_index('timestampPY',drop=True,inplace=True,verify_integrity=True) 
 

lev30cm.timestampPY = pd.to_datetime('2018-1-1') + pd.to_timedelta(lev30cm.DOY, unit='D') 
lev30cm.set_index('timestampPY',drop=True,inplace=True,verify_integrity=True)    
    


lev80cm.timestampPY = pd.to_datetime('2018-1-1') + pd.to_timedelta(lev80cm.DOY, unit='D') 
lev80cm.set_index('timestampPY',drop=True,inplace=True,verify_integrity=True)    


lev180cm.timestampPY = pd.to_datetime('2018-1-1') + pd.to_timedelta(lev180cm.DOY, unit='D') 
lev180cm.set_index('timestampPY',drop=True,inplace=True,verify_integrity=True) 

#indicate which PIC was measuring in the datafile
for df in [lev7m,lev30cm,lev80cm,lev180cm]:
    df.loc[:,'Picarro']=np.nan
    df.loc[PICbroken_2018.index[0]:PICbroken_2018.index[-1],'Picarro']='2140BROKEN'
    df.loc[PIC2120_2018.index[0]:PIC2120_2018.index[-1],'Picarro']='2120'
    df.loc[PIC2140_2018.index[0]:PIC2140_2018.index[-1],'Picarro']='2140'
#%%
exportdf(lev7m,'/Users/swa048/forServer/Vapour/2018/processed_data_combined/','EG18_level7m')
exportdf(lev30cm,'/Users/swa048/forServer/Vapour/2018/processed_data_combined/','EG18_level30cm')
exportdf(lev80cm,'/Users/swa048/forServer/Vapour/2018/processed_data_combined/','EG18_level80cm')
exportdf(lev180cm,'/Users/swa048/forServer/Vapour/2018/processed_data_combined/','EG18_level180cm')

