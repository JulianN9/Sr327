import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy import interpolate
import plotly.graph_objects as go


A = 1260*(23)**3*(10**(-3))**2
L = 0.46*10**(-3)

def rhofactor(V,R,A,L):
    # Multiplication factor:
    I = V/R
    mrho = A/(L*I)
    return mrho

def loadoldDATAlockin(path):
    cxheader = ['SignalY','SignalX','Temperature','Unknown']
    cxdata = pd.read_csv(path, skiprows=23, sep='\t', names=cxheader) 
    cxdata['SignalR'] = np.sqrt((cxdata['SignalX'])**2+(cxdata['SignalY'])**2)
    return cxdata

def loadDATAlockin(folder):
    cxpath = folder+'/'
    cxheader = ['Index','Temperature','SignalX','SignalY','Capacitancex','Capacitancey']
    cxtempdata = pd.read_csv(cxpath+'Data.csv', skiprows=1, names=cxheader)
    cxtempdata['SignalR'] = np.sqrt((cxtempdata['SignalX'])**2+(cxtempdata['SignalY'])**2) 
    return cxtempdata

def loadLARAlockin(folder):
    cxpath = folder+'/'
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

# Plotting Data
def plotdata(data, yname, save, name = ''):
    fig = plt.figure(figsize=(10,8))
    ax = fig.add_axes([0.1,0.1,0.85,0.85])
    for i in range(len(data)):
        # T = np.linspace(6,291,100)
        # Tp = np.linspace(123,291,100)
        # R = interpolate.interp1d(data[i][0]['Temperature'],data[i][2]*data[i][0][yname])
        # if i==3:
        #     print(i)
        #     ax.plot(Tp,np.gradient(np.gradient(R(Tp))),label=data[i][1],linewidth=3)
        # else:
        #     print(i)
        #     ax.plot(T,np.gradient(np.gradient(R(T))),label=data[i][1],linewidth=3)
        # if i==0:
        #     ax.plot(data[i][0]['Temperature'],0.95*data[i][2]*data[i][0][yname]+0.45,label=data[i][1],linewidth=3)        
        # elif i==2:
        #     ax.plot(data[i][0]['Temperature'],data[i][2]*data[i][0][yname]+0.45,label=data[i][1],linewidth=3)        
        # else:
        ax.plot(data[i][0]['Temperature'],data[i][2]*data[i][0][yname],label=data[i][1],linewidth=3)
    ax.set_xlabel(r'Temperature (K)',fontsize=24)
    ax.set_ylabel(r'Resistivity (m$\Omega$cm)',fontsize=24)
    ax.tick_params(axis='x', labelsize=16)
    ax.tick_params(axis='y', labelsize=16)
    # ax.set_ylabel(r'$\rho(T)$   (a.u.)',fontsize=16)
    ax.set_xlim(0,300)
    ax.legend(fontsize='16')
    if save == True:
        fig.savefig('../../Plots/Sr327/Measurements/'+name)
    plt.show()

def plotlydata(data, yname):
    fig = go.Figure()
    for i in range(len(data)):
        fig.add_trace(go.Scatter(x=data[i]['Temperature'],y=data[i][yname],name=data[i][1]))
    fig.update_layout(title='{} vs {}'.format('Temperature', 'Signal'))
    fig.update_layout(xaxis_title='{}'.format('Temperature'), yaxis_title="{}".format('Signal')) 
    fig.show()


if __name__ == "__main__":
    j_path = '../../Data/Sr327/Julian_Tom_data/'
    d_path = '../../Data/Sr327/Di_Tom_data/'

    cxcombined = combineDATAnLARA(j_path+'Version4Data/coolingdownP1',j_path+'Version4Data/coolingdownP1signal')
    cxcombinedtwo = combineDATAnLARA(j_path+'Version4Data/coolingdownP2',j_path+'Version4Data/coolingdownP2signal')
    cxcombinedthree = combineDATAnLARA(j_path+'Version3Data/heatingupP1',j_path+'Version3Data/heatingupP1signal')
    cxcombinedfour = combineDATAnLARA(j_path+'Version4Data/heatingupP1',j_path+'Version4Data/heatingupP1signal')
    cxcombinedfive = combineDATAnLARA(j_path+'Version5Data/coolingdownP1',j_path+'Version5Data/coolingdownP1signal')
    cxcombinedsix = combineDATAnLARA(j_path+'Version5Data/heatingupP1',j_path+'Version5Data/heatingupP1signal')
    cxcombinedN = combineDATAnLARA(j_path+'Version4Data/nitrogenP1',j_path+'Version4Data/nitrogenP1signal')
    ndata = pd.concat([cxcombined[(cxcombined['Temperature']>146)&(cxcombined['Temperature']<292)]+0.0000015,cxcombinedfour[(cxcombinedfour['Temperature']>=7)&(cxcombinedfour['Temperature']<=146)].iloc[::-1],cxcombined[cxcombined['Temperature']<7],cxcombinedtwo])
    # ndata = pd.concat([cxcombined[cxcombined['Temperature']>150],cxcombinedfour[(cxcombinedfour['Temperature']<=150)].iloc[::-1]])

    data = [[loadoldDATAlockin(d_path+'Sr327_Dec_05_2017_1_new.lvm'),'original c-axis',1.6],[loadDATAlockin(j_path+'Di_sample/isolation_c_cooling_down_2'),'c-axis',1.6],[loadDATAlockin(j_path+'Di_sample/isolation_d_coolingdown'),'diagonal',1.6],[loadDATAlockin(j_path+'Di_sample/isolation_x_coolingdown'),'in-plane',1.6],[ndata,'long c-axis',100*10**(-3)*rhofactor(5,1200,A,L)/83.35]]

    plotdata(data,'SignalR',True,'QualifierMeasurementFigure.svg')


    # for i in np.linspace(1.76+3*np.pi/2,1.762+3*np.pi/2,10):
    #     plotthetaskew(ndata,i)
    # plt.plot(ndata['Temperature'],ndata['SignalR'], c='blue', label='C-axis')

    # plotthetaskew(cxcombined,np.pi/4)

    # Attempting to fit a function to the new "intrinsic" value:

    # def R_c(T,A,B,C,D,E,F,G,H):  # in m-Ohm-cm (hence factor of )
    #   f = 1.0/(np.exp((T-10)/10.) + 1.)
    #   g = 1.0/(np.exp((T-40)/20.) + 1.)
    #   return (f*(A+E*T*T) + F*T*g*(B-f) + (G+H*T)*(C-g)*(D-f))

    # T = np.linspace(4.5,291,500)
    # # Tlow = np.linspace(4.5,50,100)
    # # Thigh = np.linspace(51,291,100)
    # # T = np.append(Tlow,Thigh)
    # R_i = interpolate.interp1d(ndata['Temperature'],ndata['SignalR']*100*10**(-3)*rhofactor(5,1200,A,L)/83.35)
    # R_data = R_i(T)
    # plt.scatter(T,R_data,label='interpolation', color='black')
    # guess = [1.0,1.0,1.0,1.0,0.022,0.1,7.5,0.0005*30]
    # bounds = [[0,0,0,0,0,0,5,0],[10,10,10,10,1,1,10,1]]
    # plt.plot(T,R_c(T,*guess),label='guess')

    # popt, pcov = curve_fit(R_c,T,R_data,guess,bounds=bounds)
    # print(popt)
    # plt.plot(ndata['Temperature'],R_c(ndata['Temperature'],*popt),label='fit')