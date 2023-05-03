import os
import pandas as pd
import numpy as np
from scipy import interpolate
import matplotlib.pyplot as plt
from Sr327 import loadabdata, loadcxdata2016, loadcxdata2017, plotolddataRvT

# Fitting Metallic resistivity (T^2 to T^5)
def ATsquared(T,A):
    return A*T**2

def ATn(T,A,n):
    return A*T**n

def plot2dpddataRvT(datalist, save, name = ''):
    data, labels = datalist
    T = np.linspace(6,288,100)
    for i in range(len(data)):
        R = interpolate.interp1d(data[i]['Temperature'],data[i]['Resistivity'])
        rhoc = R(T)
        plt.plot(T,np.gradient(np.gradient(rhoc)),label=labels[i])
    plt.xlim(0,300)
    plt.xlabel(r'$T$ (K)',fontsize=16)
    plt.ylabel(r'$\rho(T)$   (a.u.)',fontsize=16)
    plt.legend()
    if save == True:
        plt.savefig('../../Plots/Sr327/'+name)
    plt.show()

def dividedata(numerator, denominator):
    # XN = numerator['Temperature']; YN = numerator['Resistivity']; XD = denominator['Temperature']; YD = denominator['Resistivity']
    X = np.linspace(6,288,300)
    YD = interpolate.interp1d(denominator['Temperature'],denominator['Resistivity'])
    YN = interpolate.interp1d(numerator['Temperature'],numerator['Resistivity'])
    Y = YN(X)/YD(X)
    # plt.plot(X,YN(X))
    # plt.plot(X,YD(X))
    # plt.plot(X,Y)
    # plt.show()
    return [X, Y]

if __name__ == "__main__":
    abdata, ablabels = loadabdata()
    cxdata2016, cxlabels2016 = loadcxdata2016()
    cxdata2017, cxlabels2017 = loadcxdata2017()

    # Dividing 2016/ab:
    for number, label in enumerate(cxlabels2016):
        if label == 1.8:
            X2, Y2 = dividedata(cxdata2016[number],abdata[0])
        if label == 3.3:
            X4, Y4 = dividedata(cxdata2016[number],abdata[1])
        if label == 4.9:
            X5, Y5 = dividedata(cxdata2016[number],abdata[2])
        if label == 5.8:
            X6, Y6 = dividedata(cxdata2016[number],abdata[3])

    fig = plt.figure(figsize=(10,8))
    ax = fig.add_axes([0.1,0.1,0.85,0.85])
    ax.plot(X2,Y2,label='2Gpa')
    ax.plot(X4,Y4,label='4Gpa')
    ax.plot(X5,Y5,label='5Gpa')
    ax.plot(X6,Y6,label='6Gpa')
    ax.set_xlabel('Temperature  (K)',fontsize=16)
    ax.set_ylabel(r'$\rho_c(T)/\rho_{ab}(T)$',fontsize=16)
    
    # ax2 = fig.add_axes([0.4,0.4,0.4,0.4])
    # ax2.plot()
    # ax2.plot(X2,Y2,label='2Gpa')
    # ax2.plot(X4,Y4,label='4Gpa')
    # ax2.plot(X5,Y5,label='5Gpa')
    # ax2.plot(X6,Y6,label='6Gpa')
    # ax2.set_xlim(100,300)
    # ax2.set_ylim(0,25000)
    # ax2.set_xlabel('Temperature  (K)',fontsize=16)
    # ax2.set_ylabel(r'$\rho_c(T)/\rho_{ab}(T)$',fontsize=16)

    # ax.plot(cxdata2016[number]['Temperature'],cxdata2016[number]['Resistivity'])
    # ax.plot(abdata[0]['Temperature'],abdata[0]['Resistivity'])
    ax.legend()
    fig.savefig('../../Plots/Sr327/Rhocoverab.svg')
    plt.show()


    # plotolddataRvT(loadabdata(),False,'DiAB.svg')
    # plotolddataRvT(loadcxdata2016(),False,'DiCaxis.svg')
    # plotolddataRvT(loadcxdata2017(),False,'Di2Caxis.svg')
    # plot2dpddataRvT(loadabdata(),False,'DiABDerivative.svg')
    # plot2dpddataRvT(loadcxdata2016(),False,'DiCaxisDerivative.svg')
    # plot2dpddataRvT(loadcxdata2017(),False,'Di2CaxisDerivative.svg')