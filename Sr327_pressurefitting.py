import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import interpolate
from scipy.optimize import curve_fit
from math import floor, exp, ceil
from Sr327_poisson_simulations_burak import simulate
from Sr327 import loadcxdata2016, loadcxdata2017

def rho_c(T,P,alpha,beta,gamma,save=False,name=''):
    Nx = 24; Ny = 5
    input = 5; output = 6
    mirrorinputpin = "V["+str(Nx-input)+",1]"; mirroroutputpin = "V["+str(Nx-output)+","+str(Ny)+"]"

    df = simulate(5,6,P,alpha,beta,gamma,save=save,name=name)[0][0]

    rho = interpolate.interp1d(df['T'],2.65*(df[mirroroutputpin]-df[mirrorinputpin])/df["I"])

    Tspace = np.linspace(6,291,100)
    rhomax = max(rho(Tspace))

    return (rho(T))/rhomax

if __name__ == "__main__":
    [datalist, labellist] = loadcxdata2017()
    data = datalist[0]

    rhoc_raw = data['Resistivity']
    T_raw = data['Temperature']
    T = np.linspace(6,291,100)
    R_i = interpolate.interp1d(T_raw,rhoc_raw)
    rhodata = R_i(T)/max(R_i(T))

    guess = [labellist[0],1,1,1]
    # C_guess = rho_c(6,guess[0],guess[1],guess[2],guess[3],guess[4])-R_i(6)
    # print(C_guess)
    # guess[-1] = C_guess
    bounds = [[0,0.1,0.1,0.1],[10**(20),10,10,10]]
    params, cov = curve_fit(rho_c,T,rhodata, p0=guess,bounds=bounds)
    # fitresults = rho_c(T,params)
    # fitresults = fitresults.to_numpy()

    print(params)
    print(cov)

    fig = plt.figure(figsize=(10,8))
    ax = fig.add_axes([0.1,0.1,0.85,0.85])
    ax.plot(T,rhodata)
    ax.plot(T,rho_c(T,params[0],params[1],params[2],params[3],save=True,name='DiFitData'+str(labellist[0])+'DAT/'),ls='--')
    plt.show()

    # fig = plt.figure(figsize=(10,8))
    # ax = fig.add_axes([0.1,0.1,0.85,0.85])
    # ax.plot(T,rhodata)
    # ax.plot(T,rho_c(T,1,1,1,1,0,save=True),ls='--')
    # plt.show()