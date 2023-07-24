import os
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import interpolate
from scipy.optimize import curve_fit
from scipy.stats import chisquare
from math import floor, exp, ceil
from Sr327_poisson_simulations_burak import simulate
from Sr327_plotsimulations import RVTaxes
from Sr327 import loadcxdata2016, loadcxdata2017

def loadfitvalues(path):
    pf = pd.read_csv(path,header=[0])
    fit = pf.head(1).values.tolist()
    return fit

def rho_c(T,Pone,Ptwo,alpha,beta,gamma,Oone,Otwo,xone,xtwo,wone,wtwo,save=False,name='',conlyflag=True):
    Nx = 24; Ny = 5
    input = 5; output = 6
    mirrorinputpin = "V["+str(Nx-input)+",1]"; mirroroutputpin = "V["+str(Nx-output)+","+str(Ny)+"]"

    df = simulate(5,6,P=Pone,Ptwo=Ptwo,alpha=alpha,beta=beta,gamma=gamma,Oone=Oone,Otwo=Otwo,xone=xone,xtwo=xtwo,wone=wone,wtwo=wtwo,save=save,name=name,conlyflag=conlyflag)[0][0]

    rho = interpolate.interp1d(df['T'],2.65*(df[mirroroutputpin]-df[mirrorinputpin])/df["I"])

    Tspace = np.linspace(6,291,100)
    rhomax = max(rho(Tspace))

    return (rho(T))/rhomax

def fitpressureimplicit(data,pressure,save=True):
    #Loading previous fit data:
    path = '../../Data/Sr327_ImplicitSimulator/DiFitData'+str(pressure)+'DAT/infolist.dat'
    if Path(path).is_file() == True:
        infolist = loadfitvalues(path)[0]
    else:
        infolist = np.zeros(len(fitheader))
        infolist[4:6] = 1
        infolist[-2:] = 1
    [Pone,Ptwo,alpha,beta,gamma,Oone,Otwo,xone,xtwo,wone,wtwo] = infolist[4:]

    #Loading and formatting the Data
    rhoc_raw = data['Resistivity']
    T_raw = data['Temperature']
    T = np.linspace(6,291,100)
    # T = np.linspace(6,100,100)
    # T = np.linspace(100,291,100)
    R_i = interpolate.interp1d(T_raw,rhoc_raw)
    rhodata = R_i(T)/max(R_i(T))
    infolist[0:4] = [pressure,np.sum(R_i(T)),R_i(6),R_i(291)]

    #Fitting the data
    fittingindex = [ fitheader.index(parameter) for parameter in ["Pone","Ptwo","alpha","beta","gamma","Oone","Otwo"] ]
    guess = [infolist[j] for j in fittingindex]
    print("Guess: "+str(guess))
    lowerbounds = [10**(-10),10**(-10),-1,-1,-1,-10**3,-10**3,-5,-15,10**(-3),10**(-3)]
    upperbounds = [10**(20),10**(20),3,1,1,10**3,10**3,15,15,10**3,10**3]
    bounds = [[lowerbounds[j-4] for j in fittingindex],[upperbounds[j-4] for j in fittingindex]]
    print("Bounds: " + str(bounds))
    params, cov = curve_fit(lambda T, Pone,Ptwo,alpha,beta,gamma,Oone,Otwo: rho_c(T,Pone,Ptwo,alpha,beta,gamma,Oone,Otwo,xone,xtwo,wone,wtwo),T,rhodata, p0=guess,bounds=bounds)

    #Saving the resulting fit params
    rhofit = rho_c(T,*infolist[4:],save=True,name='DiFitData'+str(pressure)+'DAT/',conlyflag=False)
    print("Params: "+str(params))
    print("Covariance: "+str(cov))
    for count, index in enumerate(fittingindex):
        infolist[index] = params[count]

    #Plotting the resulting fit
    fig = plt.figure(figsize=(10,8))
    ax = fig.add_axes([0.1,0.1,0.85,0.85])
    ax.plot(T,rhodata,linewidth=4,label= r'Measured Data')
    ax.plot(T,rhofit,linewidth=4,ls='-.',label=r'Fitted Simulation')
    RVTaxes(ax)
    if save == True:
        fig.savefig('../../Plots/Sr327/Simulations/PressureFitting/DiFitData'+str(pressure)+'Comparison.svg')
    # plt.show()
    return infolist

if __name__ == "__main__":
    #Initalizing a global header variable:
    global fitheader 
    fitheader = ["Pressure","Scale","ScaleLow","ScaleHigh","Pone","Ptwo","alpha","beta","gamma","Oone","Otwo","xone","xtwo","wone","wtwo"]

    #Loading the data:
    [datalist, labellist] = loadcxdata2017()
    for i in range(1):
        #Initializing and Printing the Pressure Value
        data = datalist[i+5]; pressurevalue = labellist[i+5]
        print('Pressure:'+str(pressurevalue))

        #Fitting
        infolist = fitpressureimplicit(data,pressurevalue)

        #Formatting and Saving the results
        infolist = np.array(infolist).reshape(1,len(fitheader))
        df = pd.DataFrame(infolist)
        df.to_csv('../../Data/Sr327_ImplicitSimulator/DiFitData'+str(pressurevalue)+'DAT/infolist.dat', index=False, header=fitheader)

    # fig = plt.figure(figsize=(10,8))
    # ax = fig.add_axes([0.1,0.1,0.85,0.85])
    # ax.plot(T,rhodata)
    # ax.plot(T,rho_c(T,1,1,1,1,0,save=True),ls='--')
    # plt.show()