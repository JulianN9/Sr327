import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import interpolate
from scipy.optimize import curve_fit
from math import floor, exp, ceil
from Sr327_poisson_simulations_burak import simulate
from Sr327_plotsimulations import RVTaxes
from Sr327 import loadcxdata2016, loadcxdata2017

def rho_c(T,Pone,Ptwo,alpha,beta,gamma,Offset,xtwo,wtwo,save=False,name=''):
    Nx = 24; Ny = 5
    input = 5; output = 6
    mirrorinputpin = "V["+str(Nx-input)+",1]"; mirroroutputpin = "V["+str(Nx-output)+","+str(Ny)+"]"

    df = simulate(5,6,P=Pone,Ptwo=Ptwo,alpha=alpha,beta=beta,gamma=gamma,offset=Offset,xtwo=xtwo,wtwo=wtwo,save=save,name=name,conlyflag=True)[0][0]

    rho = interpolate.interp1d(df['T'],2.65*(df[mirroroutputpin]-df[mirrorinputpin])/df["I"])

    Tspace = np.linspace(6,291,100)
    rhomax = max(rho(Tspace))

    return (rho(T))/rhomax

def fitpressureimplicit(data,pressure,infolist,save=True):
    #Loading and formatting the Data
    rhoc_raw = data['Resistivity']
    T_raw = data['Temperature']
    T = np.linspace(6,291,100)
    R_i = interpolate.interp1d(T_raw,rhoc_raw)
    rhodata = R_i(T)/max(R_i(T))
    infolist.append(np.sum(R_i(T)))
    infolist.append(R_i(6))
    infolist.append(R_i(291))

    #Fitting the data
    guess = [pressure,pressure,0,0,0,0,0,1]
    bounds = [[10**(-10),10**(-10),0,-1,-1,-10**3,-15,10**(-3)],[10**(20),10**(20),3,1,1,10**3,15,10**3]]
    params, cov = curve_fit(rho_c,T,rhodata, p0=guess,bounds=bounds)

    #Saving the resulting fit params
    print(params)
    print(cov)
    for i in range(len(params)):
        infolist.append(params[i])

    #Plotting the resulting fit
    fig = plt.figure(figsize=(10,8))
    ax = fig.add_axes([0.1,0.1,0.85,0.85])
    ax.plot(T,rhodata,linewidth=4,label= r'Measured Data')
    ax.plot(T,rho_c(T,*params,save=True,name='DiFitData'+str(pressure)+'DAT/'),linewidth=4,ls='-.',label=r'Fitted Simulation')
    RVTaxes(ax)
    if save == True:
        fig.savefig('../../Plots/Sr327/Simulations/PressureFitting/DiFitData'+str(pressure)+'Comparison.svg')
    # plt.show()
    return 0

if __name__ == "__main__":
    [datalist, labellist] = loadcxdata2017()
    for i in range(7):
        #Initializing and Printing the Pressure Value
        infolist = [] 
        data = datalist[i]
        pressurevalue = labellist[i]
        print('Pressure:'+str(pressurevalue))

        #Fitting
        infolist.append(pressurevalue)
        fitpressureimplicit(data,pressurevalue,infolist)

        #Formatting and Saving the results
        header = ['Pressure','Scale','ScaleLow','ScaleHigh']
        for i in range(len(infolist)-4):
           header.append('Param '+str(i))
        infolist = np.array(infolist).reshape(1,len(infolist))
        df = pd.DataFrame(infolist)
        df.to_csv('../../Data/Sr327_ImplicitSimulator/DiFitData'+str(pressurevalue)+'DAT/infolist.dat', index=False, header=header)

    # fig = plt.figure(figsize=(10,8))
    # ax = fig.add_axes([0.1,0.1,0.85,0.85])
    # ax.plot(T,rhodata)
    # ax.plot(T,rho_c(T,1,1,1,1,0,save=True),ls='--')
    # plt.show()