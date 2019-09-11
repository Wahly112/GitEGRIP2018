#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 15 12:56:36 2019

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
tdata1 = sio.loadmat(pathm1,squeeze_me=True)
tdata2= sio.loadmat(pathm2,squeeze_me=True)
tdata3=sio.loadmat(pathm3,squeeze_me=True)


#%% create dataframes according to level
cols_f=['DOY','Valve','H2O','d18O','dD','d17O','timestampPY']
lev1=pd.DataFrame(columns=cols_f)
lev2=pd.DataFrame(columns=cols_f)
lev3=pd.DataFrame(columns=cols_f)
lev4=pd.DataFrame(columns=cols_f)


for mat in [tdata1,tdata2,tdata3]:
#level1:    
    if len(mat['mean_level1'][1,:]) <7:
        dummydf=pd.DataFrame(data=mat['mean_level1'] ,columns=['DOY','Valve','H2O','d18O','dD'])        
        lev1=lev1.append(dummydf,sort=True,ignore_index=True)        
    else:    
        dummydf=pd.DataFrame(data=mat['mean_level1'] ,columns=['DOY','Valve','timestampPY','H2O','d17O','d18O','dD'])        
        lev1=lev1.append(dummydf,sort=True,ignore_index=True)
   
#level2:    
    if len(mat['mean_level2'][1,:]) <7:
        dummydf=pd.DataFrame(data=mat['mean_level2'] ,columns=['DOY','Valve','H2O','d18O','dD'])        
        lev2=lev2.append(dummydf,sort=True,ignore_index=True)        
    else:    
        dummydf=pd.DataFrame(data=mat['mean_level2'] ,columns=['DOY','Valve','timestampPY','H2O','d17O','d18O','dD'])        
        lev2=lev2.append(dummydf,sort=True,ignore_index=True)
 #level3:    
    if len(mat['mean_level3'][1,:]) <7:
        dummydf=pd.DataFrame(data=mat['mean_level3'] ,columns=['DOY','Valve','H2O','d18O','dD'])        
        lev3=lev3.append(dummydf,sort=True,ignore_index=True)        
    else:    
        dummydf=pd.DataFrame(data=mat['mean_level3'] ,columns=['DOY','Valve','timestampPY','H2O','d17O','d18O','dD'])        
        lev3=lev3.append(dummydf,sort=True,ignore_index=True)
        
#level4:    
    if len(mat['mean_level4'][1,:]) <7:
        dummydf=pd.DataFrame(data=mat['mean_level4'] ,columns=['DOY','Valve','H2O','d18O','dD'])        
        lev4=lev4.append(dummydf,sort=True,ignore_index=True)        
    else:    
        dummydf=pd.DataFrame(data=mat['mean_level4'] ,columns=['DOY','Valve','timestampPY','H2O','d17O','d18O','dD'])        
        lev4=lev4.append(dummydf,sort=True,ignore_index=True)
      
''' reset index with datetime '''
lev1.timestampPY = pd.to_datetime('2018-1-1') + pd.to_timedelta(lev1.DOY, unit='D') - pd.Timedelta(days=1)
lev1.set_index('timestampPY',drop=True,inplace=True,verify_integrity=True) 
 

lev2.timestampPY = pd.to_datetime('2018-1-1') + pd.to_timedelta(lev2.DOY, unit='D') - pd.Timedelta(days=1)
lev2.set_index('timestampPY',drop=True,inplace=True,verify_integrity=True)    
    


lev3.timestampPY = pd.to_datetime('2018-1-1') + pd.to_timedelta(lev3.DOY, unit='D') - pd.Timedelta(days=1)
lev3.set_index('timestampPY',drop=True,inplace=True,verify_integrity=True)    


lev4.timestampPY = pd.to_datetime('2018-1-1') + pd.to_timedelta(lev4.DOY, unit='D') - pd.Timedelta(days=1)
lev4.set_index('timestampPY',drop=True,inplace=True,verify_integrity=True) 

#%%
exportdf(lev1,'/Users/swa048/forServer/Vapour/2018/processed_data_combined/','EG18_level1')
exportdf(lev2,'/Users/swa048/forServer/Vapour/2018/processed_data_combined/','EG18_level2')
exportdf(lev3,'/Users/swa048/forServer/Vapour/2018/processed_data_combined/','EG18_level3')
exportdf(lev4,'/Users/swa048/forServer/Vapour/2018/processed_data_combined/','EG18_level4')