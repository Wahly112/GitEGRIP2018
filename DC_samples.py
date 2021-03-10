# -*- coding: utf-8 -*-
"""
Created on Thu Oct 04 12:02:16 2018

@author: Sonja
"""

import os
import numpy as np
import matplotlib.pyplot as plt
import datetime as dt
import pandas as pd
plt.rcParams['text.latex.preamble']=[r"\usepackage{wasysym}"]
#path2='C:/Users/Sonja/OneDrive/Dokumente/EGRIP/EGRIP_Data/plots/'
path2='/Users/swa048/Documents/OneDrive/Dokumente/EGRIP/EGRIP_Data/plots/'
plt.close('all')


#%% import DC samples of experiment 3

path='/Users/swa048/Documents/OneDrive/Dokumente/EGRIP/EGRIP_Data/'
DCfile='EGRIP2018_DC_samples.csv'

Dc3=pd.read_csv(path+DCfile,usecols=(0,1,2,3,4,5), na_values=['NaN'],dtype={'DC_EGRIP_2018_newsampleww': object})
Dc3.rename(inplace=True, columns={
    'DC_EGRIP_2018_newsampleww': 'samplename',
    'DC_EGRIP_2018_d18O_meanww': 'd18O',
    'DC_EGRIP_2018_d18O_SDww': 'd18O_SD',
    'DC_EGRIP_2018_ddeut_meanww': 'dD',
    'DC_EGRIP_2018_ddeut_SDww': 'dD_SD',
    'DC_EGRIP_2018_dexcessww': 'dex',
})

Dc3_Ident=Dc3.samplename.str.split('_', expand=True)
Dc3_Ident.columns=['Exp','date','time','box','level']

Dc3=pd.concat([Dc3,Dc3_Ident],axis=1)
Dc3['timestamp'] = pd.to_datetime(Dc3['date'] + " " + Dc3['time'])
#Dc3.delete(106)
#Dc3.index=Dc3.timestamp

Box1=Dc3[Dc3.box=='1']
Box2=Dc3[Dc3.box=='2']
Box3=Dc3[Dc3.box=='3']
Val=Dc3[Dc3.box=='VAL']
Val.index = Val.timestamp
Val=Val.rename(columns={'level': 'DOY'})
Val.loc[:,'DOY'] = pd.Series(Val.index.dayofyear, index=Val.index)


Val.to_csv(path2+'DC_Val.txt', header=Val.columns, index=None, sep=',', mode='w')
Box1.to_csv(path2+'DC_Box1.txt', header=Box1.columns, index=None, sep=',', mode='w')
Box2.to_csv(path2+'DC_Box2.txt', header=Box2.columns, index=None, sep=',', mode='w')
Box3.to_csv(path2+'DC_Box3.txt', header=Box3.columns, index=None, sep=',', mode='w')



#%% Validation plot
dateRange = pd.date_range('07/25/2018 09:00', '07/27/2018 21:00', freq='6H')
xlabels=dateRange.format(formatter=lambda x: x.strftime('%Y-%m-%d %H:%M'))

fval, (ax1,ax2)=plt.subplots(2,sharex=True,figsize=(16,8))
ax1.plot_date(Val.timestamp,Val.d18O,color='indigo')
ax1.legend(['d18O'])
ax1.axvspan('07-25-2018 09:00', '07-26-2018 09:00', color='lightsteelblue', alpha=0.3)
ax1.axvspan('07-26-2018 09:00', '07-27-2018 09:00', color='lavender', alpha=0.3)
ax1.axvspan('07-27-2018 09:00', '07-27-2018 21:00', color='lightcyan', alpha=0.3)
ax2.plot_date(Val.timestamp,Val.dex,color='olivedrab')
ax2.axvspan('07-25-2018 09:00', '07-26-2018 09:00', color='lightsteelblue', alpha=0.3)
ax2.axvspan('07-26-2018 09:00', '07-27-2018 09:00', color='lavender', alpha=0.3)
ax2.axvspan('07-27-2018 09:00', '07-27-2018 21:00', color='lightcyan', alpha=0.3)
fval.autofmt_xdate()
ax2.legend(['d_excess'])
ax1.grid(True)
ax2.grid(True)
ax2.set_xticks(dateRange)
ax2.set_xticklabels(xlabels)
ax1.set_title('Diurnal Cycle - Validation')
fval.subplots_adjust(hspace=0.1)
fval.show()
fval.savefig(path2+"DC_Val.png")#, bbox_inches='tight')

#%% d18O plot

f18O, (ax1,ax2,ax3)=plt.subplots(3,sharex=True,figsize=(16,10))
ax1.plot_date(Box1.timestamp[Box1.level=='UL'],Box1.d18O[Box1.level=='UL'],':',color='indigo',label='UL')
ax1.plot_date(Box1.timestamp[Box1.level=='1'],Box1.d18O[Box1.level=='1'],'-.',color='darkcyan',label='1cm')
ax1.plot_date(Box1.timestamp[Box1.level=='2'],Box1.d18O[Box1.level=='2'],'-',color='darkgoldenrod',label='2cm')
ax1.plot_date(Box1.timestamp[Box1.level=='4'],Box1.d18O[Box1.level=='4'],'--',color='chocolate',label='4cm')
ax1.axvspan('07-25-2018 09:00', '07-26-2018 09:00', color='gold', alpha=0.2)
ax1.axvspan('07-26-2018 09:00', '07-27-2018 09:00', color='thistle', alpha=0.2)
ax1.axvspan('07-27-2018 09:00', '07-27-2018 21:00', color='palegreen', alpha=0.2)
ax1.set_ylabel(u'$\delta^{18}$O ‰')
ax1.set_title('Box1')
ax1.legend()
ax1.grid(True)


ax2.plot_date(Box2.timestamp[Box2.level=='UL'],Box2.d18O[Box2.level=='UL'],':',color='indigo',label='UL')
ax2.plot_date(Box2.timestamp[Box2.level=='1'],Box2.d18O[Box2.level=='1'],'-.',color='darkcyan',label='1cm')
ax2.plot_date(Box2.timestamp[Box2.level=='2'],Box2.d18O[Box2.level=='2'],'-',color='darkgoldenrod',label='2cm')
ax2.plot_date(Box2.timestamp[Box2.level=='4'],Box2.d18O[Box2.level=='4'],'--',color='chocolate',label='4cm')
ax2.axvspan('07-25-2018 09:00', '07-26-2018 09:00', color='gold', alpha=0.2)
ax2.axvspan('07-26-2018 09:00', '07-27-2018 09:00', color='thistle', alpha=0.2)
ax2.axvspan('07-27-2018 09:00', '07-27-2018 21:00', color='palegreen', alpha=0.2)
ax2.set_ylabel(u'$\delta^{18}$O ‰')
ax2.set_title('Box2')
ax2.legend()
ax2.grid(True)


ax3.plot_date(Box3.timestamp[Box3.level=='UL'],Box3.d18O[Box3.level=='UL'],':',color='indigo',label='UL')
ax3.plot_date(Box3.timestamp[Box3.level=='1'],Box3.d18O[Box3.level=='1'],'-.',color='darkcyan',label='1cm')
ax3.plot_date(Box3.timestamp[Box3.level=='2'],Box3.d18O[Box3.level=='2'],'-',color='darkgoldenrod',label='2cm')
ax3.plot_date(Box3.timestamp[Box3.level=='4'],Box3.d18O[Box3.level=='4'],'--',color='chocolate',label='4cm')
ax3.axvspan('07-25-2018 09:00', '07-26-2018 09:00', color='gold', alpha=0.2)
ax3.axvspan('07-26-2018 09:00', '07-27-2018 09:00', color='thistle', alpha=0.2)
ax3.axvspan('07-27-2018 09:00', '07-27-2018 21:00', color='palegreen', alpha=0.2)
ax3.set_title('Box3')
ax3.legend()
ax3.set_ylabel(u'$\delta^{18}$O ‰')
ax3.grid(True)
ax3.set_xticks(dateRange)
ax3.set_xticklabels(xlabels)

plt.figtext(0.25,0.9,'d18O',fontweight='bold',fontsize=18)
plt.figtext(0.25,0.08,'DOY 205',fontsize=14)
plt.figtext(0.6,0.08,'DOY 206',fontsize=14)
plt.figtext(0.8,0.08,'DOY 207',fontsize=14)
f18O.autofmt_xdate()
f18O.show()
f18O.savefig(path2+"DC_d18O.png")



#%% d-ecxess plot

fdex, (ax1,ax2,ax3)=plt.subplots(3,sharex=True,figsize=(16,10))
ax1.plot_date(Box1.timestamp[Box1.level=='UL'],Box1.dex[Box1.level=='UL'],':',color='indigo',label='UL')
ax1.plot_date(Box1.timestamp[Box1.level=='1'],Box1.dex[Box1.level=='1'],'-.',color='darkcyan',label='1cm')
ax1.plot_date(Box1.timestamp[Box1.level=='2'],Box1.dex[Box1.level=='2'],'-',color='darkgoldenrod',label='2cm')
ax1.plot_date(Box1.timestamp[Box1.level=='4'],Box1.dex[Box1.level=='4'],'--',color='chocolate',label='4cm')
ax1.axvspan('07-25-2018 09:00', '07-26-2018 09:00', color='gold', alpha=0.2)
ax1.axvspan('07-26-2018 09:00', '07-27-2018 09:00', color='thistle', alpha=0.2)
ax1.axvspan('07-27-2018 09:00', '07-27-2018 21:00', color='palegreen', alpha=0.2)
ax1.set_ylabel(u'd-excess ‰')
ax1.set_title('Box1')
ax1.legend()
ax1.grid(True)


ax2.plot_date(Box2.timestamp[Box2.level=='UL'],Box2.dex[Box2.level=='UL'],':',color='indigo',label='UL')
ax2.plot_date(Box2.timestamp[Box2.level=='1'],Box2.dex[Box2.level=='1'],'-.',color='darkcyan',label='1cm')
ax2.plot_date(Box2.timestamp[Box2.level=='2'],Box2.dex[Box2.level=='2'],'-',color='darkgoldenrod',label='2cm')
ax2.plot_date(Box2.timestamp[Box2.level=='4'],Box2.dex[Box2.level=='4'],'--',color='chocolate',label='4cm')
ax2.axvspan('07-25-2018 09:00', '07-26-2018 09:00', color='gold', alpha=0.2)
ax2.axvspan('07-26-2018 09:00', '07-27-2018 09:00', color='thistle', alpha=0.2)
ax2.axvspan('07-27-2018 09:00', '07-27-2018 21:00', color='palegreen', alpha=0.2)
ax2.set_ylabel(u'd-excess ‰')
ax2.set_title('Box2')
ax2.legend()
ax2.grid(True)


ax3.plot_date(Box3.timestamp[Box3.level=='UL'],Box3.dex[Box3.level=='UL'],':',color='indigo',label='UL')
ax3.plot_date(Box3.timestamp[Box3.level=='1'],Box3.dex[Box3.level=='1'],'-.',color='darkcyan',label='1cm')
ax3.plot_date(Box3.timestamp[Box3.level=='2'],Box3.dex[Box3.level=='2'],'-',color='darkgoldenrod',label='2cm')
ax3.plot_date(Box3.timestamp[Box3.level=='4'],Box3.dex[Box3.level=='4'],'--',color='chocolate',label='4cm')
ax3.axvspan('07-25-2018 09:00', '07-26-2018 09:00', color='gold', alpha=0.2)
ax3.axvspan('07-26-2018 09:00', '07-27-2018 09:00', color='thistle', alpha=0.2)
ax3.axvspan('07-27-2018 09:00', '07-27-2018 21:00', color='palegreen', alpha=0.2)
ax3.set_title('Box3')
ax3.legend()
ax3.set_ylabel(u'd-excess ‰')
ax3.grid(True)
ax3.set_xticks(dateRange)
ax3.set_xticklabels(xlabels)

plt.figtext(0.25,0.9,'d-excess',fontweight='bold',fontsize=18)
plt.figtext(0.25,0.08,'DOY 205',fontsize=14)
plt.figtext(0.6,0.08,'DOY 206',fontsize=14)
plt.figtext(0.8,0.08,'DOY 207',fontsize=14)
fdex.autofmt_xdate()
fdex.show()
fdex.savefig(path2+"DC_dex.png")

#%% mean d18O and d-excess
test=Dc3[Dc3.level=='UL'][Dc3.date=='20180725'][Dc3.time=='15:00']




# =============================================================================
# 
# np.savetxt(path+'ST_05cm.txt', np.concatenate((np.array(['date','DOY','d18O','std18O','D','stdD','dexcess'],dtype=str).reshape((1,7)),np.column_stack((timestamp05,values05)))), delimiter='  ',fmt='%s')  
# np.savetxt(path+'ST_1cm.txt', np.concatenate((np.array(['date','DOY','d18O','std18O','D','stdD','dexcess'],dtype=str).reshape((1,7)),np.column_stack((timestamp1,values1)))), delimiter='  ',fmt='%s')  
# np.savetxt(path+'ST_2cm.txt', np.concatenate((np.array(['date','DOY','d18O','std18O','D','stdD','dexcess'],dtype=str).reshape((1,7)),np.column_stack((timestamp2,values2)))), delimiter='  ',fmt='%s')  
# 
# =============================================================================
