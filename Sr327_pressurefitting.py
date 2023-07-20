import os
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import interpolate
from scipy.optimize import curve_fit
from math import floor, exp, ceil
from Sr327_poisson_simulations_burak import simulate
from Sr327_plotsimulations import RVTaxes
from Sr327 import loadcxdata2016, loadcxdata2017

# Best Fit Parameters:
#Header:
fitheader = ["Pressure","Scale","ScaleLow","ScaleHigh","Pone","Ptwo","Alpha","Beta","Gamma","Offset","xone","xtwo","wone","wtwo"]

def loadfitvalues(path):
    pf = pd.read_csv(path,header=[0])
    fit = pf.head(1).values.tolist()
    return fit

def rho_c(T,Pone,Ptwo,alpha,beta,gamma,Offset,xone,xtwo,wone,wtwo,save=False,name=''):
    Nx = 24; Ny = 5
    input = 5; output = 6
    mirrorinputpin = "V["+str(Nx-input)+",1]"; mirroroutputpin = "V["+str(Nx-output)+","+str(Ny)+"]"

    df = simulate(5,6,P=Pone,Ptwo=Ptwo,alpha=alpha,beta=beta,gamma=gamma,offset=Offset,xone=xone,xtwo=xtwo,wone=wone,wtwo=wtwo,save=save,name=name,conlyflag=True)[0][0]

    rho = interpolate.interp1d(df['T'],2.65*(df[mirroroutputpin]-df[mirrorinputpin])/df["I"])

    Tspace = np.linspace(6,291,100)
    rhomax = max(rho(Tspace))

    return (rho(T))/rhomax

def fitpressureimplicit(data,pressure,infolist,save=True):
    #Loading previous fit data:
    path = '../../Data/Sr327_ImplicitSimulator/DiFitData'+str(pressure)+'DAT/infolist.dat'
    if Path(path).is_file() == True:
        infolist.extend(loadfitvalues(path)[0])
    else:
        infolist = np.zeros(len(fitheader))
    [Pone,Ptwo,alpha,beta,gamma,Offset,xone,xtwo,wone,wtwo] = infolist[4:]

    #Loading and formatting the Data
    rhoc_raw = data['Resistivity']
    T_raw = data['Temperature']
    T = np.linspace(6,291,100)
    R_i = interpolate.interp1d(T_raw,rhoc_raw)
    rhodata = R_i(T)/max(R_i(T))
    # infolist.append(np.sum(R_i(T)))
    # infolist.append(R_i(6))
    # infolist.append(R_i(291))
    infolist[0:4] = [pressure,np.sum(R_i(T)),R_i(6),R_i(291)]

    #Fitting the data
    guess = [1,1]
    bounds = [[10**(-3),10**(-3)],[10**(10),10**(10)]]
    params, cov = curve_fit(lambda T, Pone, Ptwo: rho_c(T,Pone,Ptwo,alpha,beta,gamma,Offset,xone,xtwo,wone,wtwo),T,rhodata, p0=guess,bounds=bounds)

    #Saving the resulting fit params
    print(params)
    print(cov)
    # infolist.extend(fit)
    # for i in range(len(params)):
    #     infolist.insert(-2+i,params[i])
    fittingindex = [ fitheader.index(parameter) for parameter in ["Pone","Ptwo"] ]
    for count, index in enumerate(fittingindex[4:]):
        infolist[index] = params[count]

    #Plotting the resulting fit
    fig = plt.figure(figsize=(10,8))
    ax = fig.add_axes([0.1,0.1,0.85,0.85])
    ax.plot(T,rhodata,linewidth=4,label= r'Measured Data')
    ax.plot(T,rho_c(T,*infolist[4:],save=True,name='DiFitData'+str(pressure)+'DAT/'),linewidth=4,ls='-.',label=r'Fitted Simulation')
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
        # infolist.append(pressurevalue)
        # infolist[0] = pressurevalue
        fitpressureimplicit(data,pressurevalue,infolist)

        #Formatting and Saving the results
        # header = ['Pressure','Scale','ScaleLow','ScaleHigh']
        # header.extend(fitheader)
        # for i in range(len(fitheader)):
        #    header.append('Param '+str(i))
        infolist = np.array(infolist).reshape(1,len(fitheader))
        df = pd.DataFrame(infolist)
        df.to_csv('../../Data/Sr327_ImplicitSimulator/DiFitData'+str(pressurevalue)+'DAT/infolist.dat', index=False, header=fitheader)

    # fig = plt.figure(figsize=(10,8))
    # ax = fig.add_axes([0.1,0.1,0.85,0.85])
    # ax.plot(T,rhodata)
    # ax.plot(T,rho_c(T,1,1,1,1,0,save=True),ls='--')
    # plt.show()