import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from math import ceil
from pathlib import Path
import argparse
import imageio.v2 as imageio
import os
from Sr327_plotsimulations import loadsimdata, contouraxes, RVTaxes, f

if __name__ == "__main__":
    Nx = 24; Ny = 5
    input = 5; output = 6
    # p = 1
    type = 'c'
    datalist = []

    inputpin = "V["+str(input+1)+",1]"; outputpin = "V["+str(output+1)+","+str(Ny)+"]"
    mirrorinputpin = "V["+str(Nx-input)+",1]"; mirroroutputpin = "V["+str(Nx-output)+","+str(Ny)+"]"
    leftinput = "V[1,1]"; leftoutput = "V[1,"+str(Ny)+"]"
    middleinput = "V["+str(ceil(Nx/2))+",1]"; middleoutput = "V["+str(ceil(Nx/2))+","+str(Ny)+"]"
    righttinput = "V["+str(Nx)+",1]"; rightoutput = "V["+str(Nx)+","+str(Ny)+"]"

    fig = plt.figure(figsize=(10,8))
    ax = fig.add_axes([0.1,0.1,0.85,0.85])

    data, L_check = loadsimdata(type,Nx,Ny,input,output,pressure = 1,implicit = True)

    ax.scatter(data[0]["T"],0.36*data[0]["rx"],color='red',s=4.0,marker='1',label= r'$\rho_{ab}$(T)')

    for p in range(10):
        P = 1*p+1
        data, L_check = loadsimdata(type,Nx,Ny,input,output,pressure = P,implicit = True)
        ax.scatter(data[0]["T"],2.65*(data[0][mirroroutputpin]-data[0][mirrorinputpin])/data[0]["I"],s=8.0,label='C-axis, p='+str(P))

    RVTaxes(ax)
    # fig.savefig('../../Plots/Sr327/ImplicitSimulations/caxisPressureRangePFTsquared1-100.svg')
    plt.show()
    
    