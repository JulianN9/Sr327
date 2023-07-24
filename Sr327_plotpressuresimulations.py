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

    namelist = ['DiFitData1.6DAT','DiFitData2.4DAT','DiFitData3.1DAT','DiFitData3.7DAT','DiFitData4.2DAT','DiFitData4.7DAT','DiFitData5.0DAT']

    # datalist = []
    # for datatype in ['c','L']:
    #     data, L_check = loadsimdata(datatype,Nx,Ny,input,output,folder=str(Nx)+'x'+str(Ny)+'DAT',implicit = True)
    #     datalist.append(data)

    # ax.scatter(data[0]["T"],0.36*data[0]["rx"],color='red',s=4.0,marker='1',label= r'$\rho_{ab}$(T)')

    scaleheader = ['Pressure','ScaleLow','ScaleHigh','Param 0','Param 1','Param 2','Param 3','Param 4','Param 5']

    for p in range(len(namelist)-1):
        data, L_check = loadsimdata('c',Nx,Ny,input,output,implicit = True, folder=namelist[p])
        dataL, L_check = loadsimdata('L',Nx,Ny,input,output,implicit = True, folder=namelist[p])

        scalepath = '../../Data/Sr327_ImplicitSimulator/' + namelist[p] + '/infolist.dat'
        scalepd = pd.read_csv(scalepath,header=[0])
        LowTr = (data[0][mirroroutputpin][0]-data[0][mirrorinputpin][0])/data[0]["I"][0]
        HighTr = (data[0][mirroroutputpin][99]-data[0][mirrorinputpin][99])/data[0]["I"][99]
        SumTr = 0
        for i in range(100):
            SumTr = SumTr + (data[0][mirroroutputpin][i]-data[0][mirrorinputpin][i])/data[0]["I"][i]
        # scale = -(scalepd['ScaleLow'][0]-scalepd['ScaleHigh'][0])/(LowTr-HighTr)
        scale = scalepd['Scale'][0]/SumTr
        # offset = scalepd['ScaleLow'][0]-scale*LowTr
        # scale = scalepd['ScaleHigh'][0]/HighTr
        # scale = scalepd['ScaleLow'][0]/LowTr
        
        if type == 'c':
            ax.plot(data[0]["T"],scale*(data[0][mirroroutputpin]-data[0][mirrorinputpin])/data[0]["I"],label='C-axis, p='+namelist[p])
        elif type == 'L':
            ax.plot(dataL[0]["T"],-scale*(dataL[0]["V[8,1]"]-dataL[0]["V["+str(Nx-7)+",1]"])/dataL[0]["I"],label='C-axis, p='+namelist[p]) 
            ax.set_yscale("log")

    RVTaxes(ax)
    # fig.savefig('../../Plots/Sr327/Simulations/PressureFitting/'+str(type)+'axisDiDataFitted8ParamsNew.svg')
    plt.show()
    
    