#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 14 11:07:38 2019

EGRIP 2019 POWERSPECTRA of Picarro to calculate cut off frequency..

@author: swa048
"""

from loadandconvert import interpdf,deltatogas,exportdf,loadECperiod,make_patch_spines_invisible
import pandas as pd
import numpy as np
import datetime as dt
import matplotlib.pyplot as plt
from scipy import interpolate
from scipy import signal
from scipy import stats as scistats
import pwlf

#%% define function

def binaveragePS_PIC19(dfPIC,starttimes,colname='q',binnum=400,avper=30):  
    '''This function is for calculating the spectra from the half hour intervals of the Picarro continous measuremnts, interpolated frequencies 1s'''
    fmt='%m.%d %H%M'
    checkstarttime=[]
    df_fa=1/(dfPIC.index[1]-dfPIC.index[0])
    
    
#    dfPSbinned=pd.DataFrame(index=np.arange(binnum))
    dfPSbinned=pd.DataFrame()
    cutfreq=[]
    slopes=[]
    '''df = interpolated Picarro dataframe (output of interpdf())
        colname = string of column name which should be analysed for frequencies'''
    lenPIC=[]
    for idx,value in enumerate(starttimes):         
    #   interpolated and gas dataframe...
        startPIC=dfPIC.index>=starttimes[idx]
        endPIC=dfPIC.index<=starttimes[idx]+dt.timedelta(minutes=avper)
        dummyPIC=[all(tup) for tup in zip(startPIC,endPIC)]
        PIC30min=dfPIC[dummyPIC]
        lenPIC.append(len(PIC30min))
        if len(PIC30min)<df_fa*60*avper-3:   #to drop time periods with less than the needed amount of data
            continue        
        #at which time
        checkstarttime.append(PIC30min.index[0])
        #Power density spectrum
        freqs,PSden=signal.periodogram(np.array(PIC30min[colname][:1797]),fs=df_fa,detrend='linear',scaling='density')        # because of small timeshifts there can be less or more seconds in one half hour,this way, it takes only 1797 seconds
        #stats,binedges,binnumber_vec=scistats.binned_statistic(freqs[1:],PSden[1:],statistic='mean',bins=binnum)    #(x, values, statistic='mean', bins=10, range=None)  #0 frequency cut out altough mean is already subtracted because
        #dfPSbinned[start[idx].strftime(fmt)]=pd.Series(data=stats,index=dfPSbinned.index)   #
        dfPSbinned[starttimes[idx].strftime(fmt)]=pd.Series(data=PSden[1:],index=freqs[1:])
        
        x=np.log10(freqs[1:])
        y=np.log10(PSden[1:])
        #calculate cut off frequency
        my_pwlf = pwlf.PiecewiseLinFit(x, y)
        breaks = my_pwlf.fit (3)
        slope1 = my_pwlf.calc_slopes()
        
        cutfreq.append(breaks[1])
        slopes.append(slope1[0])
        
#    binmiddle=np.zeros(binnum)
#    for i in np.arange(binnum):
#        binmiddle[i]=(binedges[i]+binedges[i+1])/2   #calculate middle points of bins
#    dfPSbinned.index=binmiddle
    return dfPSbinned,checkstarttime,cutfreq,slopes
#%%
def binaveragePS_EC(dfEC,starttimes,colname='q',binnum=400,avper=30):  
    '''This function is for calculating the spectra from the half hour intervals of the Teflon tube at EGRIP19'''

    fmt='%m.%d %H%M'
    #dfPSbinned=pd.DataFrame(index=np.arange(binnum))
    dfPSbinned=pd.DataFrame()
    
    dfEC=dfEC[colname]
    
    for idx,value in enumerate(starttimes):         
    #   interpolated and gas dataframe...
        startEC=dfEC.index>=starttimes[idx]
        endEC=dfEC.index<=starttimes[idx]+dt.timedelta(minutes=avper)
        dummyEC=[all(tup) for tup in zip(startEC,endEC)]
        EC30min=dfEC[dummyEC]
        EC30min=EC30min.dropna()   #because detrend can not have nan

        if len(EC30min)<20*60*avper:   #to drop time periods with less than the needed amount of data
            continue        

        #Power density spectrum
        freqs,PSden=signal.periodogram(np.array(EC30min),fs=20,detrend='linear',scaling='density')        # because of small timeshifts there can be less or more seconds in one half hour
        #stats,binedges,binnumber_vec=scistats.binned_statistic(freqs[1:],PSden[1:],statistic='mean',bins=binnum)    #(x, values, statistic='mean', bins=10, range=None)  #0 frequency cut out altough mean is already subtracted because
        dfPSbinned[starttimes[idx].strftime(fmt)]=pd.Series(data=PSden[1:],index=freqs[1:])   #
    
#    binmiddle=np.zeros(binnum)
#    for i in np.arange(binnum):
#        binmiddle[i]=(binedges[i]+binedges[i+1])/2   #calculate middle points of bins
#    dfPSbinned.index=binmiddle
    return dfPSbinned

#%% Import the already interpolated Picarro files 

dfPICall=pd.read_csv('/Users/swa048/forServer/Vapour/2019/EGRIP19_PICinterp.txt',index_col=0,parse_dates=True,na_values=['NAN'])
# get rid of calibration data (Valve code 32)
dfPIC=dfPICall[dfPICall.Valve==4]

# transform to gas concentration instead of delta notation
deltatogas(dfPIC)   #change =17=TRUE when B for O17 is known

#import20Hz data despiked
dfIRGAq=pd.read_csv('/Users/swa048/forServer/Meteo/EC/2019/CR6/'+'EGRIP_CR6q_despiked.txt',index_col=0,parse_dates=True,na_values=['NAN'])
dfIRGA=dfcr6f

#import the 
dfLICOR=pd.read_csv('/Users/swa048/forServer/Meteo/EC/2019/CR5000/EGRIP_CR5Fluxall.txt',index_col=0,parse_dates=True,na_values=['NAN'])
dfLICOR.rename(inplace=True,columns={'H2O_mmolpm3':'q'})

#%%
#first look at data
ovview, ax1=plt.subplots()
ax1.plot(dfPIC.gas18O,'*',color='blue')
ax2=ax1.twinx()
ax2.plot(dfPIC.gasD,'s',color='brown')
ovview.legend()

#%%
loadnewPIC=0
if loadnewPIC:
    #Picarro
    PS_PIC19_q,ECstart,cutfreq_q,slopes_q=binaveragePS_PIC19(dfPIC,colname='q',avper=30)
    PS_PIC19_gas18O,ECstart,cutfreq,slopes=binaveragePS_PIC19(dfPIC,colname='gas18O',avper=30)
    PS_PIC19_gasD,ECstart,cutfreq,slopes=binaveragePS_PIC19(dfPIC,colname='gasD',avper=30)
    import csv
    pathsave='/Users/swa048/Documents/OneDrive/Dokumente/EGRIP/EGRIP2019/dataframes/'    
    PS_PIC19_q.to_csv(pathsave+'PS_PIC19_q_nb'+'.txt', header=PS_PIC19_q.columns, sep=',', mode='w',quoting=csv.QUOTE_NONNUMERIC)
    PS_PIC19_gas18O.to_csv(pathsave+'PS_PIC19_gas18O_nb'+'.txt', header=PS_PIC19_gas18O.columns, sep=',', mode='w',quoting=csv.QUOTE_NONNUMERIC)  
    PS_PIC19_gasD.to_csv(pathsave+'PS_PIC19_gasD_nb'+'.txt', header=PS_PIC19_gasD.columns, sep=',', mode='w',quoting=csv.QUOTE_NONNUMERIC)  

     
else:
    #non binned 
    file=open('/Users/swa048/Documents/OneDrive/Dokumente/EGRIP/EGRIP2019/dataframes/'+'goodtimescr6_q.txt','r')
    ECgood=file.read().split(' ')
    ECgood=pd.to_datetime(ECgood,infer_datetime_format=True)
    PS_PIC19_q=pd.read_csv('/Users/swa048/Documents/OneDrive/Dokumente/EGRIP/EGRIP2019/dataframes/PS_PIC19_q_nb.txt',index_col=0,parse_dates=True,na_values=['NAN'])
    PS_PIC19_gas18O=pd.read_csv('/Users/swa048/Documents/OneDrive/Dokumente/EGRIP/EGRIP2019/dataframes/PS_PIC19_gas18O_nb.txt',index_col=0,parse_dates=True,na_values=['NAN'])
    PS_PIC19_gasD=pd.read_csv('/Users/swa048/Documents/OneDrive/Dokumente/EGRIP/EGRIP2019/dataframes/PS_PIC19_gasD_nb.txt',index_col=0,parse_dates=True,na_values=['NAN'])
    
    
    
loadnewEC=1
if loadnewEC:
    #EC with same times
    file=open('/Users/swa048/Documents/OneDrive/Dokumente/EGRIP/EGRIP2019/dataframes/'+'goodtimescr6_q.txt','r')
    ECgoodcr6=file.read().split(' ')
    ECgoodcr6=pd.to_datetime(ECgoodcr6,infer_datetime_format=True)
    PS_IRGA19=binaveragePS_EC(dfIRGAq,colname='q',starttimes=ECgoodcr6)
    PS_IRGA19_w=binaveragePS_EC(dfIRGA,colname='Uz',starttimes=ECgoodcr6)
    file5=open('/Users/swa048/Documents/OneDrive/Dokumente/EGRIP/EGRIP2019/dataframes/'+'goodtimescr5_1.txt','r')
    ECgoodcr5=file5.read().split(' ')
    ECgoodcr5=pd.to_datetime(ECgoodcr5,infer_datetime_format=True)
    PS_LICOR19=binaveragePS_EC(dfLICOR,starttimes=ECgoodcr5)
    #PS_IRGA19_detrend=binaveragePS_EC(dfIRGA,ECstart[100:200])
    #save dataframes
    import csv
    pathsave='/Users/swa048/Documents/OneDrive/Dokumente/EGRIP/EGRIP2019/dataframes/'
    PS_IRGA19.to_csv(pathsave+'PS_IRGA19_nb'+'.txt', header=PS_IRGA19.columns, sep=',', mode='w',quoting=csv.QUOTE_NONNUMERIC) 
    PS_IRGA19_w.to_csv(pathsave+'PS_IRGA19_w_nb'+'.txt', header=PS_IRGA19_w.columns, sep=',', mode='w',quoting=csv.QUOTE_NONNUMERIC) 
    PS_LICOR19.to_csv(pathsave+'PS_LICOR19_nb'+'.txt', header=PS_LICOR19.columns, sep=',', mode='w',quoting=csv.QUOTE_NONNUMERIC) 
    #PS_IRGA19_detrend.to_csv(pathsave+'PS_IRGA19_detrend'+'.txt', header=PS_IRGA19_detrend.columns, sep=',', mode='w',quoting=csv.QUOTE_NONNUMERIC) 
else:
    PS_IRGA19=pd.read_csv('/Users/swa048/Documents/OneDrive/Dokumente/EGRIP/EGRIP2019/dataframes/PS_IRGA19_nb.txt',index_col=0,parse_dates=True,na_values=['NAN'])
    #PSIRGA_binned=pd.read_csv('/Users/swa048/Documents/OneDrive/Dokumente/EGRIP/EGRIP2019/dataframes/PS_IRGA19.txt',index_col=0,parse_dates=True,na_values=['NAN'])
    PS_IRGA19_w=pd.read_csv('/Users/swa048/Documents/OneDrive/Dokumente/EGRIP/EGRIP2019/dataframes/PS_IRGA19_w_nb.txt',index_col=0,parse_dates=True,na_values=['NAN'])
    PS_LICOR19=pd.read_csv('/Users/swa048/Documents/OneDrive/Dokumente/EGRIP/EGRIP2019/dataframes/PS_LICOR19_nb.txt',index_col=0,parse_dates=True,na_values=['NAN'])
    #PS_IRGA19_detrend=pd.read_csv('/Users/swa048/Documents/OneDrive/Dokumente/EGRIP/EGRIP2019/dataframes/PS_IRGA19_detrend.txt',index_col=0,parse_dates=True,na_values=['NAN']) #a reduced number of ogive computations,linealary detrended instead of constat detrend
    #save dataframes


#%%  calculate mean
PS_PIC19_q['qmean']=PS_PIC19_q.mean(axis=1)
PS_PIC19_gas18O['gas18Omean']=PS_PIC19_gas18O.mean(axis=1)
PS_PIC19_gasD['gasDmean']=PS_PIC19_gasD.mean(axis=1)
PS_IRGA19['qmean']=PS_IRGA19.mean(axis=1)
#PS_IRGA19_detrend['qmean']=PS_IRGA19_detrend.mean(axis=1)
PS_IRGA19_w['wmean']=PS_IRGA19_w.mean(axis=1)
PS_LICOR19['qmean']=PS_LICOR19.mean(axis=1)
#%%  plot the different Averaged Power spectra
means,ax1=plt.subplots(figsize=(12,6))
ax1.loglog(PS_PIC19_q.index,PS_PIC19_q.qmean,color='darkcyan',label='q mean constant')
ax1.loglog(PS_PIC19_q.index,PS_PIC19_q[PS_PIC19_q.columns[200]],color='darkcyan',alpha=0.5)


#h=ax1.get_ylim()
#logdiff=np.log10(h[1])-np.log10(h[0])
#ax4=ax1.twinx()
ax1.loglog(PS_IRGA19.index,PS_IRGA19.qmean,color='darkorchid',label='q mean IRGA')
#ax4.spines["left"].set_position(("axes", 1.1))
#make_patch_spines_invisible(ax4)
#ax4.spines["left"].set_visible(True)
#ax2=ax1.twinx()
ax1.loglog(PS_PIC19_gas18O.index,PS_PIC19_gas18O.gas18Omean,color='gray',label='d18O')
ax1.loglog(PS_PIC19_gas18O.index,PS_PIC19_gas18O[PS_PIC19_gas18O.columns[200]],color='gray',alpha=0.5)
#h2=ax2.get_ylim()
#newhigh=10**(np.log10(h2[0])+logdiff)
#ax2.set_ylim(h2[0],newhigh)
#ax2.tick_params(axis='y', which='minor')
#ax2.yaxis.set_minor_formatter(plt.LogFormatter())
#
#ax3=ax1.twinx()
ax1.loglog(PS_PIC19_gasD.index,PS_PIC19_gasD.gasDmean,color='olivedrab',label='dD')
ax1.loglog(PS_PIC19_gasD.index,PS_PIC19_gasD[PS_PIC19_gas18O.columns[200]],color='olivedrab',alpha=0.5)
#ax3.spines["right"].set_position(("axes", 1.07))
#make_patch_spines_invisible(ax3)
#ax3.spines["right"].set_visible(True)
#h3=ax3.get_ylim()
#newhigh3=10**(np.log10(h3[0])+logdiff)
#ax3.set_ylim(h3[0],newhigh3)
ax1.grid(True)
means.legend()
#%%  try to find the slope of the mean, first one slope, then piecewise

#find slope of q (the most important) then for 18gas etc-- 
slopesq,ax1=plt.subplots(figsize=(12,6))
ax1.loglog(PS_PIC19_q.index,PS_PIC19_q.qmean,color='darkcyan',label='q mean')

# linear fit model for loglog plot
x=np.log10(PS_PIC19_q.index)
y=np.log10(PS_PIC19_q.qmean)

m,c=np.polyfit(x,y,1)
yfit=10**(m*x+c)
#ax1.loglog(PS_PIC19_q.index,yfit,':',label='lin plot')
print(m)

# piecewise linear fit with pwlf
my_pwlf = pwlf.PiecewiseLinFit(x, y)
breaksq = my_pwlf.fit (3)
slopes = my_pwlf.calc_slopes()
print('the first slope is'+ str(slopes[0]))
print(breaksq)
ynew=10**(my_pwlf.predict(x))
ax1.set_title('Picarro water mixing ratio signal')
ax1.set_xlabel('frequencies in Hz')
ax1.set_ylabel('PS in mixing ratio/s^2')
ax1.loglog(PS_PIC19_q.index,ynew,'r',label='piecewise plot')
ax1.text(0.4, 0.85, f'slope of first segment= {np.round(slopes[0],2)}', horizontalalignment='center',verticalalignment='center', transform=ax1.transAxes)
ax1.text(0.6, 0.65, f'the first break is at period= {np.round(1/(10**breaksq[1]))}s', horizontalalignment='center',verticalalignment='center', transform=ax1.transAxes)
ax1.loglog(PS_PIC19_q.index,1e-5*PS_PIC19_q.index**(-5/3),color='cadetblue',LineStyle=':',label='-5/3')
ax1.legend()
slopesq.savefig('/Users/swa048/Documents/OneDrive/Dokumente/EC/Isoflux/Isoflux_plots/cutfreq/'+'EGRIP19_q'+'.pdf')
#%%

# plot for gas 18O
slopes18,ax1=plt.subplots(figsize=(12,6))
ax1.loglog(PS_PIC19_gas18O.index,PS_PIC19_gas18O.gas18Omean,color='gray',label='gas18O mean')

# linear fit model for loglog plot
x=np.log10(PS_PIC19_gas18O.index)
y=np.log10(PS_PIC19_gas18O.gas18Omean)

m,c=np.polyfit(x,y,1)
yfit=10**(m*x+c)
#ax1.loglog(PS_PIC19_gas18O.index,yfit,':',label='lin plot')
print(m)

# piecewise linear fit with pwlf
my_pwlf = pwlf.PiecewiseLinFit(x, y)
breaks18 = my_pwlf.fit (3)
slopes = my_pwlf.calc_slopes()
print('the first gas 18O slope is'+ str(slopes[0]))
print(breaks18)

ynew=10**(my_pwlf.predict(x))
ax1.set_title('Picarro d18O gas signal')
ax1.set_xlabel('frequencies in Hz')
ax1.set_ylabel('PS in mixing ratio/s^2')

ax1.loglog(PS_PIC19_gas18O.index,ynew,'r',label='piecewise plot')
ax1.text(0.4, 0.85, f'slope of first segment= {np.round(slopes[0],2) }', horizontalalignment='center',verticalalignment='center', transform=ax1.transAxes)
ax1.text(0.6, 0.65, f'the first break is at period= {np.round(1/(10**breaks18[1])) }s', horizontalalignment='center',verticalalignment='center', transform=ax1.transAxes)

ax1.loglog(PS_PIC19_gas18O.index,1e-10*PS_PIC19_gas18O.index**(-5/3),color='cadetblue',LineStyle=':',label='-5/3')
ax1.legend()
slopes18.savefig('/Users/swa048/Documents/OneDrive/Dokumente/EC/Isoflux/Isoflux_plots/cutfreq/'+'EGRIP19_18O'+'.pdf')

#%%
# plot for gas Deuterium
slopesD,ax1=plt.subplots(figsize=(12,6))
ax1.loglog(PS_PIC19_gasD.index,PS_PIC19_gasD.gasDmean,color='olivedrab',label='gasD mean')

# linear fit model for loglog plot
x=np.log10(PS_PIC19_gasD.index)
y=np.log10(PS_PIC19_gasD.gasDmean)

m,c=np.polyfit(x,y,1)
yfit=10**(m*x+c)
#ax1.loglog(PS_PIC19_gasD.index,yfit,':',label='lin plot')
print(m)

# piecewise linear fit with pwlf
my_pwlf = pwlf.PiecewiseLinFit(x, y)
breaksD = my_pwlf.fit (3)
slopes = my_pwlf.calc_slopes()
print('the first gas D slope is'+ str(slopes[0]))
print(breaksD)

ynew=10**(my_pwlf.predict(x))
ax1.set_title('Picarro dD gas signal')
ax1.set_xlabel('frequencies in Hz')
ax1.set_ylabel('PS in mixing ratio/s^2')

ax1.loglog(PS_PIC19_gasD.index,ynew,'r',label='piecewise plot')
ax1.text(0.4, 0.85, f'slope of first segment= {np.round(slopes[0],2)}', horizontalalignment='center',verticalalignment='center', transform=ax1.transAxes)
ax1.text(0.6, 0.65, f'the first break is at period= {np.round(1/(10**breaksD[1]))}s', horizontalalignment='center',verticalalignment='center', transform=ax1.transAxes)
ax1.loglog(PS_PIC19_gasD.index,1e-13*PS_PIC19_gasD.index**(-5/3),color='cadetblue',LineStyle=':',label='-5/3')
ax1.legend()
slopesD.savefig('/Users/swa048/Documents/OneDrive/Dokumente/EC/Isoflux/Isoflux_plots/cutfreq/'+'EGRIP19_D'+'.pdf')
#%%
# plot for humidity from Irgason in g/m3
slopesIRGA,ax1=plt.subplots(figsize=(12,6))
ax1.loglog(PS_IRGA19.index,PS_IRGA19.qmean,color='maroon',label='Irga hum mean despiked')


# linear fit model for loglog plot
xb=np.log10(PS_IRGA19.index)
yb=np.log10(PS_IRGA19.qmean)

my_pwlf_b = pwlf.PiecewiseLinFit(xb, yb)
breaksIRGA_b = my_pwlf_b.fit (2)
slopes_b = my_pwlf_b.calc_slopes()
print('the first hum slope for binned freqs is'+ str(slopes_b[0]))
print(breaksIRGA_b)
ynew_b=10**(my_pwlf_b.predict(xb))

ax1.loglog(PS_IRGA19.index,ynew_b,'pink',label='piecewise plot')
ax1.text(0.4, 0.85, f'slope of first segment= {np.round(slopes_b[0],2)}', horizontalalignment='center',verticalalignment='center', transform=ax1.transAxes)
ax1.text(0.6, 0.65, f'the first break is at period= {np.round(1/(10**breaksIRGA_b[1]))}s', horizontalalignment='center',verticalalignment='center', transform=ax1.transAxes)
ax1.loglog(PS_IRGA19.index,1e-5*PS_IRGA19.index**(-5/3),color='cadetblue',LineStyle=':',label='-5/3')
ax1.set_title('Irgason humidity signal')
ax1.set_xlabel('frequencies in Hz')
ax1.set_ylabel('PS in (g/m^3)/s^2')
ax1.legend()
slopesIRGA.savefig('/Users/swa048/Documents/OneDrive/Dokumente/EC/Isoflux/Isoflux_plots/cutfreq/'+'EGRIP19_IRGA'+'.pdf')
#%% Plot for LICOR data
slopesLICOR,ax1=plt.subplots(figsize=(12,6))
ax1.loglog(PS_LICOR19.index,PS_LICOR19.qmean,color='maroon',label='LICOR hum mean')


# linear fit model for loglog plot
x=np.log10(PS_LICOR19.index)
y=np.log10(PS_LICOR19.qmean)


m,c=np.polyfit(x,y,1)
yfit=10**(m*x+c)
ax1.loglog(PS_LICOR19.index,yfit,':',color='orange',label='lin plot')
print(m)

# piecewise linear fit with pwlf
my_pwlf = pwlf.PiecewiseLinFit(x, y)
breaksLICOR = my_pwlf.fit (3)
slopes = my_pwlf.calc_slopes()
print('the first hum slope is'+ str(slopes[0]))
print(breaksLICOR)
ynew=10**(my_pwlf.predict(x))


ax1.loglog(PS_LICOR19.index,ynew,'r',label='piecewise plot_binned')
ax1.text(0.4, 0.85, f'slope of first segment= {np.round(slopes[0],2)}', horizontalalignment='center',verticalalignment='center', transform=ax1.transAxes)
ax1.text(0.6, 0.65, f'the first break is at period= {np.round(1/(10**breaksLICOR[1]))}s', horizontalalignment='center',verticalalignment='center', transform=ax1.transAxes)
ax1.loglog(PS_LICOR19.index,PS_LICOR19.index**(-5/3),color='cadetblue',LineStyle=':',label='-5/3')
ax1.set_title('LICOR humidity signal')
ax1.set_xlabel('frequencies in Hz')
ax1.set_ylabel('PS in (mol/m^3)/s^2')
ax1.legend()
slopesLICOR.savefig('/Users/swa048/Documents/OneDrive/Dokumente/EC/Isoflux/Isoflux_plots/cutfreq/'+'EGRIP19_LICOR'+'.pdf')
#%% IRGA_w plot
slopes_w,ax1=plt.subplots(figsize=(12,6))
ax1.loglog(PS_IRGA19_w.index,PS_IRGA19_w.wmean,color='maroon',label='vertical wind')


# linear fit model for loglog plot
x=np.log10(PS_IRGA19_w.index)
y=np.log10(PS_IRGA19_w.wmean)


m,c=np.polyfit(x,y,1)
yfit=10**(m*x+c)
ax1.loglog(PS_IRGA19_w.index,yfit,':',color='orange',label='lin plot')
print(m)

# piecewise linear fit with pwlf
my_pwlf = pwlf.PiecewiseLinFit(x, y)
breaks_w = my_pwlf.fit (3)
slopes = my_pwlf.calc_slopes()
print('the first hum slope is'+ str(slopes[0]))
print(breaks_w)
ynew=10**(my_pwlf.predict(x))


ax1.loglog(PS_IRGA19_w.index,ynew,'r',label='piecewise plot_binned')
ax1.text(0.4, 0.85, f'slope of first segment= {np.round(slopes[0],2)}', horizontalalignment='center',verticalalignment='center', transform=ax1.transAxes)
ax1.text(0.6, 0.65, f'the first break is at period= {np.round(1/(10**breaks_w[1]))}s', horizontalalignment='center',verticalalignment='center', transform=ax1.transAxes)
ax1.loglog(PS_IRGA19_w.index,1e-2*PS_IRGA19_w.index**(-5/3),color='cadetblue',LineStyle=':',label='-5/3')
ax1.set_title('vertical wind signal')
ax1.set_xlabel('frequencies in Hz')
ax1.set_ylabel('PS in (m/s)/s^2')
ax1.legend()
slopes_w.savefig('/Users/swa048/Documents/OneDrive/Dokumente/EC/Isoflux/Isoflux_plots/cutfreq/'+'EGRIP19_w'+'.pdf')
#%% plot the individual PS
highfreqmean=PS_IRGA19[PS_IRGA19.columns[:-1]][5:10].mean()   
highfreqmedian=PS_IRGA19[PS_IRGA19.columns[:-1]][5:10].median()

singleps,ax1=plt.subplots()

ax1.loglog(PS_IRGA19.index,PS_IRGA19[str(highfreqmean.idxmax())],color='red') 

for col2 in PS_IRGA19.columns[:-1][highfreqmean>8*10**(-5)]:  #plot the half hours that are corrupt in light orange
    ax1.loglog(PS_IRGA19.index,PS_IRGA19[col2],color='orange',alpha=0.2)  
for col in PS_IRGA19.columns[:-1][highfreqmean<8*10**(-5)]:  #plot the good PS 
    ax1.loglog(PS_IRGA19.index,PS_IRGA19[col],color='black',alpha=0.1)

#%% plot the individual PS LICOR
highfreqmeanLICOR=PS_LICOR19[PS_LICOR19.columns[:-1]][5:10].mean()   
highfreqmedianLICOR=PS_LICOR19[PS_LICOR19.columns[:-1]][5:10].median()

singleps,ax1=plt.subplots()

ax1.loglog(PS_LICOR19.index,PS_LICOR19[str(highfreqmeanLICOR.idxmax())],color='red') 

for col2 in PS_LICOR19.columns[:-1][highfreqmeanLICOR>8*10**(-5)]:  #plot the half hours that are corrupt in light orange
    ax1.loglog(PS_LICOR19.index,PS_LICOR19[col2],color='orange',alpha=0.2)  
for col in PS_LICOR19.columns[:-1][highfreqmeanLICOR<8*10**(-5)]:  #plot the good PS 
    ax1.loglog(PS_LICOR19.index,PS_LICOR19[col],color='black',alpha=0.1)
#%% histogram of high freq means and median to filter our which PS are not usable
# the criterion is, that the mean value in the high frequencies is too high i.e. those half hour PS influence the overall mean PS and leads to the wrong shape


histo=plt.figure()

#plt.hist(np.log10(highfreqmean),100)
plt.hist(highfreqmean.dropna(),1000)
plt.xlim(0,0.001)
plt.title('mean')

histo2=plt.figure()

#plt.hist(np.log10(highfreqmedian),100)
plt.hist(highfreqmedian.dropna(),1000)
plt.xlim(0,0.001)
plt.title('median')


#%% Non Usable half hours of Irga
idx_right=PS_IRGA19.columns[:-1][highfreqmean<8*10**(-5)]
right30min= pd.to_datetime('2019.'+idx_right,format='%Y.%m.%d %H%M')
pararight=pd.DataFrame(index=right30min,columns=['Ts','wd','ws','varq'])

idx_false=PS_IRGA19.columns[:-1][highfreqmean>8*10**(-5)]
false30min= pd.to_datetime('2019.'+idx_false,format='%Y.%m.%d %H%M')
parafalse=pd.DataFrame(index=false30min,columns=['Ts','wd','ws','varq'])

computenew=0
whatswrong,(axwrong,axright,axboth)=plt.subplots(3,1,sharey=True,figsize=(16,10))
xaxis=np.arange(36001)


for hh in false30min:
    dfhh=dfIRGA.q[hh:hh+dt.timedelta(minutes=30)]
    if computenew:
        parafalse.Ts[hh]=(dfIRGA.T_SONIC[hh:hh+dt.timedelta(minutes=30)]).mean()
        parafalse.varq[hh]=dfhh.var()
        parafalse.wd[hh]=np.round((np.arctan2(-dfIRGA.Uy[hh:hh+dt.timedelta(minutes=30)],dfIRGA.Ux[hh:hh+dt.timedelta(minutes=30)])*180/np.pi+245).mean(),1)
        parafalse.ws[hh]=np.round((np.sqrt(dfIRGA.Ux[hh:hh+dt.timedelta(minutes=30)]**2+dfIRGA.Uy[hh:hh+dt.timedelta(minutes=30)]**2)).mean(),1)
        parafalse=parafalse.apply(pd.to_numeric)
    axwrong.plot(xaxis,dfhh)
    axwrong.set_title('all 170 time series of bad PS')
    axboth.plot(xaxis,dfhh,color='r',alpha=0.4)

#whatswrong.legend()

a=0
for hhyes in right30min:
    try:
        dfhhyes=dfIRGA.q[hhyes:hhyes+dt.timedelta(minutes=30)]
    except KeyError:    #in case some Piucarro times are not resplved by the Irgaosn example 20.06.2019 Irgason keine Daten
        a=a+1
        continue
    if computenew:
        pararight.Ts[hhyes]=(dfIRGA.T_SONIC[hhyes:hhyes+dt.timedelta(minutes=30)]).mean()
        pararight.varq[hhyes]=dfhhyes.var()
        pararight.wd[hhyes]=np.round((np.arctan2(-dfIRGA.Uy[hhyes:hhyes+dt.timedelta(minutes=30)],dfIRGA.Ux[hhyes:hhyes+dt.timedelta(minutes=30)])*180/np.pi+245).mean(),1)
        pararight.ws[hhyes]=np.round((np.sqrt(dfIRGA.Ux[hhyes:hhyes+dt.timedelta(minutes=30)]**2+dfIRGA.Uy[hhyes:hhyes+dt.timedelta(minutes=30)]**2)).mean(),1)
        pararight=pararight.apply(pd.to_numeric)
    axright.plot(xaxis,dfhhyes,label=str(hhyes))
    axright.set_title('all 641 time series of good PS')
    axboth.plot(xaxis,dfhhyes,color='k',alpha=0.2)
#thatsright.legend()    


axwrong.set_ylim(0,25)
axright.set_ylim(0,25)
axboth.set_ylim(0,25)
axwrong.set_ylabel('Irga humidity in g/m3')
axright.set_ylabel('Irga humidity in g/m3')
axboth.set_ylabel('Irga humidity in g/m3')

whatswrong.savefig('/Users/swa048/Documents/OneDrive/Dokumente/EC/Isoflux/Isoflux_plots/cutfreq/'+'EGRIP19_whatswrong'+'.png')
#%% look at parameters to identify causes for bad timeseries
import windrose
import matplotlib.cm as cm

parahist,((ax1,ax2,ax3),(ax4,ax5,ax6),)=plt.subplots(2,3,figsize=(16,20))
#ax1=windspeed
ax1.hist(pararight.dropna().ws, 100,color='green',label='no spikes' )
ax1.hist(parafalse.ws , 100,color='red',label='spikes' )
ax1.legend()
ax1.set_xlabel('m/s')
ax1.set_title('wind speed')
#ax2=varq
ax2.hist(pararight.dropna().varq, 100,color='green',label='no spikes' )
ax2.hist(parafalse.varq , 100,color='red',label='spikes' )
ax2.legend()
ax2.set_title('variance q')
#ax3=varq
ax3.hist(pararight.dropna().Ts, 100,color='green',label='no spikes' )
ax3.hist(parafalse.Ts , 100,color='red',label='spikes' )
ax3.legend()
ax3.set_title('sonic Temp')
ax3.set_xlabel('degree D')
#ax4=windrose
ax4.remove()
ax4=parahist.add_subplot(2,3,4,projection='windrose')
ax4.bar(parafalse.wd-245, parafalse.ws, normed=True, opening=1, edgecolor='white',cmap=cm.winter)
ax4.set_title('spikes')
ax5.remove()
ax5=parahist.add_subplot(2,3,5,projection='windrose')
ax5.bar(pararight.wd-245, pararight.ws, normed=True, opening=1, edgecolor='white',cmap=cm.winter)
ax5.set_title('no spikes')

parahist.savefig('/Users/swa048/Documents/OneDrive/Dokumente/EC/Isoflux/Isoflux_plots/cutfreq/'+'EGRIP19_spikesparams'+'.pdf')

#%% look at both mean PS and hopefully see a nice -5/3 slope for the good PSem
badPS=PS_IRGA19[PS_IRGA19.columns[:-1][highfreqmean>8*10**(-5)]]
badPSqmean=badPS.mean(axis=1)

goodPS=PS_IRGA19[PS_IRGA19.columns[:-1][highfreqmean<8*10**(-5)]]
goodPSqmean=goodPS.mean(axis=1)

compPS=plt.figure()
plt.loglog(badPSqmean,label='badPS')
plt.loglog(goodPSqmean,label='goodPS')
compPS.legend()
#%%   not really useful

##normalize good PS
#def nearest(items, pivot):
#    return min(items, key=lambda x: abs(x - pivot))
#normfreq=nearest(goodPS.index,10**-2)
##normidx=goodPS.index.get_loc(normfreq)
#
#goodPSnorm=pd.DataFrame(index=goodPS.index)
#
#for col in goodPS.columns:
#    goodPSnorm[col]=goodPS[col]/goodPS.loc[normfreq,col]
#
#goodPSnormmean=goodPSnorm.mean(axis=1)
#    
#PSnorm=plt.figure()
#for col in goodPSnorm:
#    plt.loglog(goodPSnorm.index,goodPSnorm[col],color='k',alpha=0.2)
#plt.loglog(goodPSnormmean,color='r')
#plt.loglog(goodPSnorm.index,goodPSnorm.index**(-5/3),':')
##PSnorm.legend()



#%% fitting a line through the good PS

goodPS,ax1=plt.subplots(figsize=(12,6))
ax1.loglog(goodPSqmean,color='gray',label='good PS mean')

# linear fit model for loglog plot
x=np.log10(goodPSqmean.index)
y=np.log10(goodPSqmean)

m,c=np.polyfit(x,y,1)
yfit=10**(m*x+c)
#ax1.loglog(PS_PIC19_gas18O.index,yfit,':',label='lin plot')
print(m)

xc=[x[0]]
yc=[y.iloc[0]]
# piecewise linear fit with pwlf
my_pwlf = pwlf.PiecewiseLinFit(x, y)
breaks18 = my_pwlf.fit (2,y_c=yc,x_c=xc)
slopes = my_pwlf.calc_slopes()
print('the first slope for good PS is'+ str(slopes[0]))
print(breaks18)

ynew=10**(my_pwlf.predict(x))
ax1.set_title('Irgason humidity signal')
ax1.set_xlabel('frequencies in Hz')
ax1.set_ylabel('PS')

ax1.loglog(goodPSqmean.index,ynew,'r',label='piecewise plot')
ax1.text(0.4, 0.85, f'slope of first segment= {np.round(slopes[0],2) }', horizontalalignment='center',verticalalignment='center', transform=ax1.transAxes)
ax1.text(0.6, 0.65, f'the first break is at Hz= {10**breaks18[1]}', horizontalalignment='center',verticalalignment='center', transform=ax1.transAxes)

ax1.loglog(goodPSqmean.index,1e-5*goodPSqmean.index**(-5/3),color='cadetblue',LineStyle=':',label='-5/3')
ax1.legend()
#goodPS.savefig('/Users/swa048/Documents/OneDrive/Dokumente/EC/Isoflux/Isoflux_plots/cutfreq/'+'EGRIP19_IRGAgood'+'.pdf')

