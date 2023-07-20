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

# Best Fit Parameters:
#Header:
fitheader = ["Pone","Ptwo","Beta","Gamma","Offset","xone","xtwo","wone","wtwo"]
#1.6
fit16 = [2.79969501975309,502747.9338725326,0.6451196931186023,0.10142682831433505,-0.12863793191504097,9.811356676399319,-1.1576332126515747,0.7469557326731626]
fit24 = [8.430456673694383,1877869.2330952075,0.9024214033254677,0.17013024945521896,-0.09514846054340166,11.640980396656033,-1.8846811780899666,0.7617640249824599]
fit31 = [68.44719489068251,1.6976851540898308,0.2988192167638081,0.13847485314892985,-0.03166898083181844,2.4333065107838814,-2.426596301510242,0.8967675783057795
]
fit37 = [11150075.261443684,1.212530988478126,0.1699064163343642,0.15852928370656194,0.001535463431130828,7.812992031222445,-2.631494023919922,0.7485219226053272]
fit42 = [42731939.76194198,3.6814019775230338,0.14093747228534798,0.3789907141609256,0.006486515723140473,18.006447259208205,-5.457881350102343,0.7353040915272021]
fit47 = [576.7280054343599,3.1663065622885593,0.09983839416882699,0.18511973924800892,0.000938895314810774,0.28599009182791535,-2.6953956000045696,0.8187686816131754]
fit50 = [5282.799452247552,3.1502891399845065,0.1405046777807265,0.1903764110378832,0.019147027356810987,0.6540446772444144,-2.705100973409597,0.8101864240921575]
fits = [fitheader,fit16,fit24,fit31,fit37,fit42,fit47,fit50]

def Rho(Pone,Ptwo,alpha,beta,gamma,Offset,xtwo,wtwo):
    def rho_c(T,xone,wone,save=False,name=''):
        Nx = 24; Ny = 5
        input = 5; output = 6
        mirrorinputpin = "V["+str(Nx-input)+",1]"; mirroroutputpin = "V["+str(Nx-output)+","+str(Ny)+"]"

        df = simulate(5,6,P=Pone,Ptwo=Ptwo,alpha=alpha,beta=beta,gamma=gamma,offset=Offset,xone=xone,xtwo=xtwo,wone=wone,wtwo=wtwo,save=save,name=name,conlyflag=True)[0][0]

        rho = interpolate.interp1d(df['T'],2.65*(df[mirroroutputpin]-df[mirrorinputpin])/df["I"])

        Tspace = np.linspace(6,291,100)
        rhomax = max(rho(Tspace))

        return (rho(T))/rhomax
    return rho_c

def fitpressureimplicit(data,pressure,infolist,count,save=True):
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
    guess = [0,1]
    bounds = [[-5,10**(-3)],[15,10**3]]
    params, cov = curve_fit(Rho(*fits[count]),T,rhodata, p0=guess,bounds=bounds)

    #Saving the resulting fit params
    print(params)
    print(cov)
    for i in range(len(params)):
        infolist.append(params[i])

    #Plotting the resulting fit
    fig = plt.figure(figsize=(10,8))
    ax = fig.add_axes([0.1,0.1,0.85,0.85])
    ax.plot(T,rhodata,linewidth=4,label= r'Measured Data')
    ax.plot(T,Rho(*fits[count])(T,*params,save=True,name='DiFitData'+str(pressure)+'DAT/'),linewidth=4,ls='-.',label=r'Fitted Simulation')
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
        fitpressureimplicit(data,pressurevalue,infolist,i+1)

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