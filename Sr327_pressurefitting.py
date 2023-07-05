import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import interpolate
from math import floor, exp, ceil
from Sr327_poisson_simulations_burak import simulate

def rho_c(T,P):
    Nx = 24; Ny = 5
    input = 5; output = 6
    inputpin = "V["+str(input+1)+",1]"; outputpin = "V["+str(output+1)+","+str(Ny)+"]"
    mirrorinputpin = "V["+str(Nx-input)+",1]"; mirroroutputpin = "V["+str(Nx-output)+","+str(Ny)+"]"
    leftinput = "V[1,1]"; leftoutput = "V[1,"+str(Ny)+"]"
    middleinput = "V["+str(ceil(Nx/2))+",1]"; middleoutput = "V["+str(ceil(Nx/2))+","+str(Ny)+"]"
    righttinput = "V["+str(Nx)+",1]"; rightoutput = "V["+str(Nx)+","+str(Ny)+"]"

    df = simulate(5,6,P,save=False)[0][0]

    # print(type(df))
    # fig = plt.figure(figsize=(10,8))
    # ax = fig.add_axes([0.1,0.1,0.85,0.85])
    # ax.scatter(df["T"],2.65*(df[mirroroutputpin]-df[mirrorinputpin])/df["I"],s=8.0,label='C-axis, p='+str(P))
    # plt.show()

    # I = [ i for i in range(100)]
    # I.reverse()
    # rho = 2.65*(df[mirroroutputpin][T]-df[mirrorinputpin][T])/df["I"][T]
    rho = interpolate.interp1d(df['T'],2.65*(df[mirroroutputpin]-df[mirrorinputpin])/df["I"])

    return rho(T)

if __name__ == "__main__":
    T = np.linspace(2,300,100)
    # print(T)
    plt.plot(T,rho_c(T,1))
    plt.show()    