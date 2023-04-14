import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy import interpolate


A = 1260*(23)**3*(10**(-3))**2
L = 0.46*10**(-3)

def rhofactor(V,R,A,L):
    # Multiplication factor:
    I = V/R
    mrho = A/(L*I)
    return mrho

def loadDATAlockin(folder):
    cxpath = '../../Data/'+folder+'/'
    cxheader = ['Index','Temperature','SignalX','SignalY','Capacitancex','Capacitancey']
    cxtempdata = pd.read_csv(cxpath+'Data.csv', skiprows=1, names=cxheader) 
    return cxtempdata

def loadLARAlockin(folder):
    cxpath = '../../Data/'+folder+'/'
    cxheader = ['Index','SignalX','SignalY','Time']
    cxdata = pd.read_csv(cxpath+'data_cropped.csv', skiprows=1, names=cxheader) 
    return cxdata

def combineDATAnLARA(foldertemp,foldersignal):
    cxtemp = loadDATAlockin(foldertemp)
    cxdata = loadLARAlockin(foldersignal)
    header = ['Index','Temperature','SignalX','SignalY','SignalR','Theta']
    cxcombined = pd.DataFrame(columns=header)
    cxcombined['Index'] = cxtemp['Index']
    cxcombined['Temperature'] = cxtemp['Temperature']
    cxcombined['SignalX'] = cxdata['SignalX']
    cxcombined['SignalY'] = cxdata['SignalY']
    cxcombined['SignalR'] = np.sqrt((cxdata['SignalX'])**2+(cxdata['SignalY'])**2)
    cxcombined['Theta'] = np.arctan2(cxdata['SignalY'],cxdata['SignalX'])
    return cxcombined

def plotthetaskew(dataframe,theta):
    plt.plot(dataframe['Temperature'],np.cos(theta)*dataframe['SignalX']+np.sin(theta)*dataframe['SignalY'],label=str(theta))

SR327path = 'Sr327/Julian_Tom_data/'

cxcombined = combineDATAnLARA(SR327path+'Version4Data/coolingdownP1',SR327path+'Version4Data/coolingdownP1signal')
cxcombinedtwo = combineDATAnLARA(SR327path+'Version4Data/coolingdownP2',SR327path+'Version4Data/coolingdownP2signal')
cxcombinedthree = combineDATAnLARA(SR327path+'Version3Data/heatingupP1',SR327path+'Version3Data/heatingupP1signal')
cxcombinedfour = combineDATAnLARA(SR327path+'Version4Data/heatingupP1',SR327path+'Version4Data/heatingupP1signal')
cxcombinedfive = combineDATAnLARA(SR327path+'Version5Data/coolingdownP1',SR327path+'Version5Data/coolingdownP1signal')
cxcombinedsix = combineDATAnLARA(SR327path+'Version5Data/heatingupP1',SR327path+'Version5Data/heatingupP1signal')
cxcombinedN = combineDATAnLARA(SR327path+'Version4Data/nitrogenP1',SR327path+'Version4Data/nitrogenP1signal')
ndata = pd.concat([cxcombined[(cxcombined['Temperature']>146)&(cxcombined['Temperature']<292)]+0.0000015,cxcombinedfour[(cxcombinedfour['Temperature']>=7)&(cxcombinedfour['Temperature']<=146)].iloc[::-1],cxcombined[cxcombined['Temperature']<7],cxcombinedtwo])
# ndata = pd.concat([cxcombined[cxcombined['Temperature']>150],cxcombinedfour[(cxcombinedfour['Temperature']<=150)].iloc[::-1]])
# plt.plot(cxcombinedthree['Temperature'],cxcombinedthree['SignalR'], c='black',label='Heating Up 1')
# plt.plot(cxcombinedfour['Temperature'],cxcombinedfour['SignalR'], c='red',label='Heating Up 2',ls='--')
# plt.plot(cxcombined['Temperature'],cxcombined['SignalR'], c='blue', label='Cooling Down',ls='--')
# plt.plot(cxcombinedtwo['Temperature'],cxcombinedtwo['SignalR'], c='blue',ls='--')
# plt.plot(cxcombinedfive['Temperature'],cxcombinedfive['SignalR']*100*10**(-3)*rhofactor(5,1200,A,L), c='orange',ls='--',label = 'No amp Cooling Down')
# plt.plot(cxcombinedsix['Temperature'],cxcombinedsix['SignalR']*100*10**(-3)*rhofactor(5,1200,A,L), c='black',ls='--',label='No amp Heating Up')
plt.plot(ndata['Temperature'],ndata['SignalR']*100*10**(-3)*rhofactor(5,1200,A,L)/83.35, c='blue', label='C-axis')
# plt.plot(ndata['Temperature'],ndata['SignalR']/83.35, c='blue', label='C-axis')
# plt.plot(cxcombinedN['Temperature'],cxcombinedN['SignalR']*100*10**(-3)*rhofactor(5,1200,A,L), c='grey', label='Nitrogen',ls='--')


# for i in np.linspace(1.76+3*np.pi/2,1.762+3*np.pi/2,10):
#     plotthetaskew(ndata,i)
# plt.plot(ndata['Temperature'],ndata['SignalR'], c='blue', label='C-axis')

# plotthetaskew(cxcombined,np.pi/4)

# Attempting to fit a function to the new "intrinsic" value:

def R_c(T,B,D,E,F,G,H):  # in m-Ohm-cm (hence factor of )
  f = 1.0/(np.exp((T-10)/B) + 1.)
  g = 1.0/(np.exp((T-40)/D) + 1.)
  return (f*(1.0+E*T*T) + F*T*g*(1.-f) + (G+H*T)*(1.-g)*(1.-f))


# T = np.linspace(4.5,291,100)
Tlow = np.linspace(4.5,50,100)
Thigh = np.linspace(51,291,100)
T = np.append(Tlow,Thigh)
R_i = interpolate.interp1d(ndata['Temperature'],ndata['SignalR']*100*10**(-3)*rhofactor(5,1200,A,L)/83.35)
R_data = R_i(T)
plt.scatter(T,R_data,label='interpolation', color='black')
guess = [10.0,20.0,0.022,0.1,7.5,0.0005*30]
bounds = [[1,1,0,0,7,0],[20,80,1,1,8,1]]
plt.plot(T,R_c(T,*guess),label='guess')


popt, pcov = curve_fit(R_c,T,R_data,guess,bounds=bounds)
print(popt)
plt.plot(ndata['Temperature'],R_c(ndata['Temperature'],*popt),label='fit')

# plt.plot(cxcombinedtwo['Temperature'],cxcombinedtwo['SignalY']+0.005)
# cxcombinedb = combinedata('coolingdown1V1_P2','coolingdown_signal1V1_P2')
# plt.plot(cxcombinedb['Temperature'],cxcombinedb['SignalX'])
# cxcombinedc = combinedata('coolingdown1V1_P3','coolingdown_signal1V1_P3')
# plt.plot(cxcombinedc['Temperature'],cxcombinedc['SignalX'])

plt.xlabel(r'Temperature (K)',fontsize=16)
plt.ylabel(r'Resistivity (m$\Omega$cm)',fontsize=16)
# plt.ylabel(r'$\rho(T)$   (a.u.)',fontsize=16)
plt.xlim(0,300)
plt.legend()
# plt.savefig('../../Plots/Sr327/caxis_rho_final_V1.svg')
plt.show()