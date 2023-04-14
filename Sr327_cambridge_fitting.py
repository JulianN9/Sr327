import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Fitting Metallic resistivity (T^2 to T^5)
def ATsquared(T,A):
    return A*T**2

def ATn(T,A,n):
    return A*T**n

# ab Data from Di from cambridge
def loadabdata():
    abpath = '../../data/Sr327_ab_plane_data_Cambridge/uptoroomtemp/sample2/filtered/'
    abheader = ['Temperature','Resistivity']
    abdata = []; ablabels = []
    for file in os.listdir(os.getcwd()+'/'+abpath):
        abdata.append(pd.read_csv(abpath+file, sep=' ', names=abheader))
        ablabels.append(file)
    for data in abdata:
        data[abheader[0]] = data[abheader[0]]/1000
    ablabels = [2,4,5,6,7]
    ablabels, abdata = zip(*sorted(zip(ablabels,abdata)))
    return [abdata, ablabels]

# c-axis Data from Di, 2016 is first sample, 2017 is second sample
def loadcxdata2016():
    cxpath = '../../data/Di_Tom_data/filtered/'
    cxheader = ['Unknown','Resistivity','Temperature']
    cxdata = []; cxlabels = []
    for file in os.listdir(os.getcwd()+'/'+cxpath):
        if '2016' in file:
            cxdata.append(pd.read_csv(cxpath+file, skiprows=23, sep='\t', names=cxheader))
            cxlabels.append(file)
            if 'Sep' in file:
                cxdata[-1]['Resistivity'] = cxdata[-1]['Resistivity']/10
    for data in cxdata:
        data['Resistivity'] = np.abs(data['Resistivity'])
    cxlabels = [4.9,3.3,1.8,2.8,5.8]
    cxlabels, cxdata = zip(*sorted(zip(cxlabels,cxdata)))
    return [cxdata, cxlabels]

def loadcxdata2017():
    cxpath = '../../data/Di_Tom_data/filtered/'
    cxheader = ['Unknown','Resistivity','Temperature','AlsoUnknown']
    cxdata = []; cxlabels = []
    for file in os.listdir(os.getcwd()+'/'+cxpath):
        if '2017' in file:
            cxdata.append(pd.read_csv(cxpath+file, skiprows=23, sep='\t', names=cxheader))
            cxlabels.append(file)
    for data in cxdata:
        data['Resistivity'] = np.abs(data['Resistivity'])
    cxlabels = [3.7,1.6,5.0,4.7,3.1,2.4,4.2]
    cxlabels, cxdata = zip(*sorted(zip(cxlabels,cxdata)))
    return [cxdata, cxlabels]

def plotpddataRvT(datalist, save, name = ''):
    data, labels = datalist
    for i in range(len(data)):
        plt.plot(data[i]['Temperature'],data[i]['Resistivity'],label=labels[i])
    # line = np.linspace(0,300)
    # plt.plot(line, 10**(-5)*line**2)
    # plt.ylim(0,0.05)
    plt.xlim(0,300)
    plt.xlabel(r'$T$ (K)',fontsize=16)
    plt.ylabel(r'$\rho(T)$   (a.u.)',fontsize=16)
    plt.legend()
    if save == True:
        plt.savefig(name)
    plt.show()

def dividedata(numerator, denominator):
    XN = numerator['Temperature']; YN = numerator['Resistivity']; XD = denominator['Temperature']; YD = denominator['Resistivity']
    X = []; Y = []
    if len(XD)>len(XN):
        X = XD
        YN = np.interp(X,XN,YN,left=6,right=260)
        YD = np.interp(X,XD,YD,left=6,right=260)
    else:
        X = XN
        YD = np.interp(X,XD,YD,left=6,right=260)
        YN = np.interp(X,XN,YN,left=6,right=260)
    Y = YN/YD
    # plt.plot(X,Y)
    # plt.show()
    return [X, Y]

abdata, ablabels = loadabdata()
cxdata2016, cxlabels2016 = loadcxdata2016()
cxdata2017, cxlabels2017 = loadcxdata2017()

#Dividing 2016/ab:
for number, label in enumerate(cxlabels2016):
    if label == 1.8:
        X2, Y2 = dividedata(cxdata2016[number],abdata[0])
        # plt.plot(X2,Y2/(350000/4))
        # plt.plot(cxdata2016[number]['Temperature'],cxdata2016[number]['Resistivity'])
        # plt.plot(abdata[0]['Temperature'],abdata[0]['Resistivity'])
        # plt.show()
    if label == 3.3:
        X4, Y4 = dividedata(cxdata2016[number],abdata[1])
    if label == 4.9:
        X5, Y5 = dividedata(cxdata2016[number],abdata[2])
    if label == 5.8:
        X6, Y6 = dividedata(cxdata2016[number],abdata[3])

# plotpddataRvT(loadabdata(),False)
# plotpddataRvT(loadcxdata2016(),False)
plotpddataRvT(loadcxdata2017(),False)