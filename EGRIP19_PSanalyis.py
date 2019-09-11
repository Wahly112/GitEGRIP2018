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

def binaveragePS_PIC19(dfPIC,colname='q',binnum=400,avper=30):  
    '''This function is for calculating the spectra from the half hour intervals of the Teflon tube at KOHNEN'''
    
    df_dt=(dfPIC.index[1]-dfPIC.index[0]).total_seconds()
    df_fa=round(1/df_dt,1)
    
    checkstarttime=[]
    
    
    timedt=np.zeros(len(dfPIC.index)-1)   #one shorter than original
    for i in np.arange(1,len(dfPIC.index)):     #start from 1 until end
        dummy=dfPIC.index[i]-dfPIC.index[i-1]
        timedt[i-1]=dummy.total_seconds()
    
    
    PIC_dtlog=timedt!=df_dt                     #find out all time when Valve Code 4/0 started by scanning for longer timeintervals than the interpolated df_dt interval
    PIC_dtlog=np.insert(PIC_dtlog,0,False)   #the array which to add to, the index, the value you want to add    --> because boolean area is one too short
    fmt='%m.%d %H%M'
    start=dfPIC.index[PIC_dtlog]            #the array with the onset of Valve 4/0
    
    
#    dfPSbinned=pd.DataFrame(index=np.arange(binnum))
    dfPSbinned=pd.DataFrame()
    cutfreq=[]
    slopes=[]
    '''df = interpolated Picarro dataframe (output of interpdf())
        colname = string of column name which should be analysed for frequencies'''
    lenPIC=[]
    for idx,value in enumerate(start):         
    #   interpolated and gas dataframe...
        startPIC=dfPIC.index>=start[idx]
        endPIC=dfPIC.index<=start[idx]+dt.timedelta(minutes=avper)
        dummyPIC=[all(tup) for tup in zip(startPIC,endPIC)]
        PIC30min=dfPIC[dummyPIC]
        lenPIC.append(len(PIC30min))
        if len(PIC30min)<df_fa*60*avper-3:   #to drop time periods with less than the needed amount of data
            continue        
        #at which time
        checkstarttime.append(PIC30min.index[0])
        #Power density spectrum
        freqs,PSden=signal.periodogram(np.array(PIC30min[colname][:1797]),fs=df_fa,detrend='constant',scaling='density')        # because of small timeshifts there can be less or more seconds in one half hour
        #stats,binedges,binnumber_vec=scistats.binned_statistic(freqs[1:],PSden[1:],statistic='mean',bins=binnum)    #(x, values, statistic='mean', bins=10, range=None)  #0 frequency cut out altough mean is already subtracted because
        #dfPSbinned[start[idx].strftime(fmt)]=pd.Series(data=stats,index=dfPSbinned.index)   #
        dfPSbinned[start[idx].strftime(fmt)]=pd.Series(data=PSden[1:],index=freqs[1:])
        
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
    
    for idx,value in enumerate(starttimes):         
    #   interpolated and gas dataframe...
        startEC=dfEC.index>=starttimes[idx]
        endEC=dfEC.index<=starttimes[idx]+dt.timedelta(minutes=avper)
        dummyEC=[all(tup) for tup in zip(startEC,endEC)]
        EC30min=dfEC[dummyEC]

        if len(EC30min)<20*60*avper:   #to drop time periods with less than the needed amount of data
            continue        

        #Power density spectrum
        freqs,PSden=signal.periodogram(np.array(EC30min[colname]),fs=20,detrend='linear',scaling='density')        # because of small timeshifts there can be less or more seconds in one half hour
        #stats,binedges,binnumber_vec=scistats.binned_statistic(freqs[1:],PSden[1:],statistic='mean',bins=binnum)    #(x, values, statistic='mean', bins=10, range=None)  #0 frequency cut out altough mean is already subtracted because
        dfPSbinned[starttimes[idx].strftime(fmt)]=pd.Series(data=PSden[1:],index=freqs[1:])   #
    
#    binmiddle=np.zeros(binnum)
#    for i in np.arange(binnum):
#        binmiddle[i]=(binedges[i]+binedges[i+1])/2   #calculate middle points of bins
#    dfPSbinned.index=binmiddle
    return dfPSbinned
#%%



#%% Import the already interpolated Picarro files 

dfPICall=pd.read_csv('/Users/swa048/forServer/Vapour/2019/EGRIP19_PICinterp.txt',index_col=0,parse_dates=True,na_values=['NAN'])

#import20Hz data for comparison=
dfIRGA=pd.read_csv('/Users/swa048/forServer/Meteo/EC/2019/CR6/EGRIP_CR6Fluxall.txt',index_col=0,parse_dates=True,na_values=['NAN'])
dfIRGA.rename(inplace=True,columns={'H2O_density':'q'})
# get rid of calibration data (Valve code 32)
dfPIC=dfPICall[dfPICall.Valve==4]

# transform to gas concentration instead of delta notation
deltatogas(dfPIC)   #change =17=TRUE when B for O17 is known
#%%
#first look at data
ovview, ax1=plt.subplots()
ax1.plot(dfPIC.gas18O,'*',color='blue')
ax2=ax1.twinx()
ax2.plot(dfPIC.gasD,'s',color='brown')
ovview.legend()

#%%
loadnew=0
if loadnew:
    #Picarro
    PS_PIC19_q,ECstart,cutfreq_q,slopes_q=binaveragePS_PIC19(dfPIC,colname='q',avper=30)
    PS_PIC19_gas18O,ECstart,cutfreq,slopes=binaveragePS_PIC19(dfPIC,colname='gas18O',avper=30)
    PS_PIC19_gasD,ECstart,cutfreq,slopes=binaveragePS_PIC19(dfPIC,colname='gasD',avper=30)
    
    #EC with same times
    PS_IRGA19=binaveragePS_EC(dfIRGA,ECstart)
    PS_IRGA19_detrend=binaveragePS_EC(dfIRGA,ECstart[100:200])
else:
    #non binned 
    PS_PIC19_q=pd.read_csv('/Users/swa048/Documents/OneDrive/Dokumente/EGRIP/EGRIP2019/dataframes/dataframesPS_PIC19_q_nb.txt',index_col=0,parse_dates=True,na_values=['NAN'])
    PS_PIC19_gas18O=pd.read_csv('/Users/swa048/Documents/OneDrive/Dokumente/EGRIP/EGRIP2019/dataframes/dataframesPS_PIC19_gas18O_nb.txt',index_col=0,parse_dates=True,na_values=['NAN'])
    PS_PIC19_gasD=pd.read_csv('/Users/swa048/Documents/OneDrive/Dokumente/EGRIP/EGRIP2019/dataframes/dataframesPS_PIC19_gasD_nb.txt',index_col=0,parse_dates=True,na_values=['NAN'])
    PS_IRGA19=pd.read_csv('/Users/swa048/Documents/OneDrive/Dokumente/EGRIP/EGRIP2019/dataframes/dataframesPS_IRGA19_nb.txt',index_col=0,parse_dates=True,na_values=['NAN'])
    PSIRGA_binned=pd.read_csv('/Users/swa048/Documents/OneDrive/Dokumente/EGRIP/EGRIP2019/dataframes/dataframesPS_IRGA19.txt',index_col=0,parse_dates=True,na_values=['NAN'])
    PS_IRGA19_detrend=binaveragePS_EC(dfIRGA,ECstart[100:200]) #a reduced number of ogive computations,linealary detrended instead of constat detrend
#%%  calculate mean
PS_PIC19_q['qmean']=PS_PIC19_q.mean(axis=1)
PS_PIC19_gas18O['gas18Omean']=PS_PIC19_gas18O.mean(axis=1)
PS_PIC19_gasD['gasDmean']=PS_PIC19_gasD.mean(axis=1)
PS_IRGA19['qmean']=PS_IRGA19.mean(axis=1)
PS_IRGA19_detrend['qmean']=PS_IRGA19_detrend.mean(axis=1)
#%%  plot the different Averaged Power spectra
means,ax1=plt.subplots(figsize=(12,6))
ax1.loglog(PS_PIC19_q.index,PS_PIC19_q.qmean,color='darkcyan',label='q mean constant')
ax1.loglog(PS_PIC19_q.index,PS_PIC19_q[PS_PIC19_q.columns[200]],color='darkcyan',alpha=0.5)


ax1.loglog(PS_IRGA19_detrend.index,PS_IRGA19_detrend.qmean,color='darkcyan',label='q mean linear')

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
slopesD,ax1=plt.subplots(figsize=(12,6))
ax1.loglog(PS_IRGA19_detrend.index,PS_IRGA19_detrend.qmean,color='blueviolet',label='Irga hum mean linear')
ax1.loglog(PS_IRGA19.index,PS_IRGA19.qmean,color='maroon',label='Irga hum mean const')


# linear fit model for loglog plot
x=np.log10(PS_IRGA19_detrend.index)
y=np.log10(PS_IRGA19_detrend.qmean)

xb=np.log10(PS_IRGA19.index)
yb=np.log10(PS_IRGA19.qmean)

m,c=np.polyfit(x,y,1)
yfit=10**(m*x+c)
ax1.loglog(PS_IRGA19_detrend.index,yfit,':',color='orange',label='lin plot')
print(m)

# piecewise linear fit with pwlf
my_pwlf = pwlf.PiecewiseLinFit(x, y)
breaksIRGA = my_pwlf.fit (3)
slopes = my_pwlf.calc_slopes()
print('the first hum slope is'+ str(slopes[0]))
print(breaksIRGA)
ynew=10**(my_pwlf.predict(x))


my_pwlf_b = pwlf.PiecewiseLinFit(xb, yb)
breaksIRGA_b = my_pwlf_b.fit (2)
slopes_b = my_pwlf_b.calc_slopes()
print('the first hum slope for binned freqs is'+ str(slopes_b[0]))
print(breaksIRGA_b)
ynew_b=10**(my_pwlf_b.predict(xb))

ax1.loglog(dfIRGA.index,ynew_b,'pink',label='piecewise plot_binned')
ax1.loglog(PS_IRGA19_detrend.index,ynew,'r',label='piecewise plot_binned')
ax1.text(0.4, 0.85, f'slope of first segment= {np.round(slopes[0],2)}', horizontalalignment='center',verticalalignment='center', transform=ax1.transAxes)
ax1.text(0.6, 0.65, f'the first break is at period= {np.round(1/(10**breaksIRGA[1]))}s', horizontalalignment='center',verticalalignment='center', transform=ax1.transAxes)
ax1.loglog(PS_IRGA19_detrend.index,1e-5*PS_IRGA19_detrend.index**(-5/3),color='cadetblue',LineStyle=':',label='-5/3')
ax1.set_title('Irgason humidity signal')
ax1.set_xlabel('frequencies in Hz')
ax1.set_ylabel('PS in (g/m^3)/s^2')
ax1.legend()
#slopesD.savefig('/Users/swa048/Documents/OneDrive/Dokumente/EC/Isoflux/Isoflux_plots/cutfreq/'+'EGRIP19_IRGA'+'.pdf')
#%% plot the individual PS
highfreqmean=PS_IRGA19[PS_IRGA19.columns[:-1]][5:10].mean()   
highfreqmedian=PS_IRGA19[PS_IRGA19.columns[:-1]][5:10].median()

singleps,ax1=plt.subplots()

ax1.loglog(PS_IRGA19.index,PS_IRGA19[str(highfreqmean.idxmax())],color='red') 

for col in PS_IRGA19.columns[:-1]:  #plot all different half hours in light grey
    ax1.loglog(PS_IRGA19.index,PS_IRGA19[col],color='black',alpha=0.1)
for col2 in PS_IRGA19.columns[:-1][highfreqmean>8*10**(-4)]:  #plot the half hours that are corrupt in light orange
    ax1.loglog(PS_IRGA19.index,PS_IRGA19[col2],color='orange',alpha=0.2)  
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
#%% save the dataframes
if loadnew:
    pathsave='/Users/swa048/Documents/OneDrive/Dokumente/EGRIP/EGRIP2019/dataframes/'
    exportdf(PS_PIC19_q,pathsave,'PS_PIC19_q_nb')    
    exportdf(PS_PIC19_gas18O,pathsave,'PS_PIC19_gas18O_nb')  
    exportdf(PS_PIC19_gasD,pathsave,'PS_PIC19_gasD_nb')  
    exportdf(PS_IRGA19,pathsave,'PS_IRGA19_nb')
#%% Non Usable half hours of Irga
idx_false=PS_IRGA19.columns[:-1][highfreqmean>8*10**(-4)]
false30min= pd.to_datetime('2019.'+idx_false,format='%Y.%m.%d %H%M')

whatswrong=plt.figure()
xaxis=np.arange(36001)

for hh in false30min:
    dfhh=dfIRGA.q[hh:hh+dt.timedelta(minutes=30)]
    plt.plot(xaxis,dfhh,label=str(hh))

whatswrong.legend()

idx_right=PS_IRGA19.columns[:100][highfreqmean[:100]<8*10**(-4)]
right30min= pd.to_datetime('2019.'+idx_right,format='%Y.%m.%d %H%M')

thatsright=plt.figure()
xaxis=np.arange(36001)

for hhyes in right30min:
    dfhhyes=dfIRGA.q[hhyes:hhyes+dt.timedelta(minutes=30)]
    plt.plot(xaxis,dfhhyes,label=str(hhyes))

thatsright.legend()
#%% look at both mean PS and hopefully see a nice -5/3 slope for the good PSem
badPS=PS_IRGA19[PS_IRGA19.columns[:-1][highfreqmean>8*10**(-4)]]
badPS['qmean']=badPS.mean(axis=1)

goodPS=PS_IRGA19[PS_IRGA19.columns[:-1][highfreqmean<8*10**(-4)]]
goodPS['qmean']=goodPS.mean(axis=1)

compPS=plt.figure()
plt.plot(badPS.index,badPS.qmean,label='badPS')
plt.plot(goodPS.index,goodPS.qmean,label='goodPS')
compPS.legend()