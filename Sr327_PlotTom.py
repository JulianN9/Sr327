from asyncore import write
from cProfile import label
from xml.dom.minidom import Element
from importlib_metadata import files
from requests import head
from sklearn import mixture
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.mlab
import pandas as pd
import scipy
import plotly.express as px
from scipy import signal
from scipy.signal import argrelextrema
from pylab import *
from scipy.optimize import leastsq
import plotly.graph_objects as go
import csv

#Reading Data
path = '../Data/Julian_Sr327_Tom_data_2022/'
header = [['SignalY','SignalX','Temperature','Unknown'],['Index','Temperature','SignalX','SignalY','Capacitancex','Capacitancey']]
# files = [['Sr327_Dec_05_2017_1_new.lvm',0,'Original c-axis'],['isolation_c_cooling_down_2/Data.csv',1,'New Iso c-axis'],['isolation_d_coolingdown/Data.csv',1,'New Iso Diagonal'],['isolation_x_coolingdown/Data.csv',1,'New Iso In-Plane']]
# files = [['Sr327_Dec_05_2017_1_new.lvm',0,'Original c-axis'],['isolation_c_cooling_down_2/Data.csv',1,'New New C long']]
files = [[path+'Sr327_Dec_05_2017_1_new.lvm',0,'Original c-axis'],[path+'Disample1_c/Data.csv',1,'New c-axis'],[path+'Disample1_diagonal_coolingdown/Data.csv',1,'New Diagonal'],[path+'Disample1_inplane_coolingdown/Data.csv',1,'New In-Plane']]
heating = False
if heating == True:
    files.append([path+'Disample1_c_heatingup/Data.csv',1,'New Heating Up'])
    files.append([path+'Disample1_diagonal_heatingup/Data.csv',1,'Diagonal Heating Up'])
    files.append([path+'Disample1_inplane_heatingup/Data.csv',1,'In-Plane Heating Up'])
data = []
for file in files:
    if file[1] == 0:
        data.append(pd.read_csv(file[0], skiprows=23, sep='\t', names=header[0])) 
    if file[1] == 1:
        data.append(pd.read_csv(file[0], skiprows=1, names=header[1]))

#print(data[0]['Temperature'].to_numpy())

# templist = data[3]['Temperature'].to_numpy()
# xlist = data[3]['SignalX'].to_numpy()
# ylist = data[3]['SignalY'].to_numpy()
# savenumb = 0

# for i in range(0,66061-1):
#     if np.sqrt(xlist[i]**2+ylist[i]**2) > 1.5:
#         if templist[i] > savenumb:
#             print(templist[i])
#             savenumb = templist[i]

# Plotting Data
def plotdata(data, save, name = ''):
    for i in range(len(data)):
        plt.plot(data[i][header[0][2]],np.sqrt(data[i][header[0][1]]**2+data[i][header[0][0]]**2),label=files[i][2])
    plt.xlabel(r'$T$ (K)',fontsize=16)
    plt.ylabel(r'$\rho(T)$   (a.u.)',fontsize=16)
    plt.legend()
    if save == True:
        plt.savefig(name)
    plt.show()

def plotlydata(data):
    fig = go.Figure()
    for i in range(len(data)):
        fig.add_trace(go.Scatter(x=data[i][header[0][2]],y=np.sqrt(data[i][header[0][1]]**2+data[i][header[0][0]]**2),name=files[i][2]))
    fig.update_layout(title='{} vs {}'.format(header[0][2], 'Signal'))
    fig.update_layout(xaxis_title='{}'.format(header[0][2]), yaxis_title="{}".format('Signal')) 
    fig.show()

plotdata(data, False, 'DiFig_Clean.svg')
#plotlydata(data)