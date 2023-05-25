import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy import interpolate
import plotly.graph_objects as go
import Sr327
from Sr327_plotsimulations import loadsimdata, RVTaxes
from math import ceil

Sr327.sr327globalsetup()

data = Sr327.loadqualifierdata()
Nx = 24
Ny = 5
input = 5
output = 6
inputpin = "V["+str(input+1)+",1]"; outputpin = "V["+str(output+1)+","+str(Ny)+"]"
mirrorinputpin = "V["+str(Nx-input)+",1]"; mirroroutputpin = "V["+str(Nx-output)+","+str(Ny)+"]"
leftinput = "V[1,1]"; leftoutput = "V[1,"+str(Ny)+"]"
middleinput = "V["+str(ceil(Nx/2))+",1]"; middleoutput = "V["+str(ceil(Nx/2))+","+str(Ny)+"]"
righttinput = "V["+str(Nx)+",1]"; rightoutput = "V["+str(Nx)+","+str(Ny)+"]"
simdata, L_check = loadsimdata('m',Nx,Ny,5,6,0)


# ax = fig.add_axes([0.1,0.1,0.85,0.85])
fig = plt.figure(figsize=(16,7))
(ax1, ax2) = fig.subplots(1,2)
Sr327.plotdata(fig,ax1,data,'SignalR',False,'QualifierMeasurementFigure.svg',False,True,original=False)
ax2.plot(simdata[0]["T"],2.65*(simdata[0][mirroroutputpin]-simdata[0][mirrorinputpin])/simdata[0]["I"],color='C1',marker='^',label='c-axis')
ax2.plot(simdata[1]["T"],2.65*(simdata[1][outputpin]-simdata[1][mirrorinputpin])/simdata[1]["I"],color='green',marker='^',label='diagonal')
ax2.plot(simdata[2]["T"],(Ny/(Nx-output-input))*2.65*(simdata[2][mirroroutputpin]-simdata[2][outputpin])/simdata[2]["I"],color='red',marker='^',label='in-plane')
#ax.plot(simdata[2]["T"],2.65*(simdata[2][mirrorinputpin]-simdata[2][inputpin])/simdata[2]["I"],color='red',marker='^')
if L_check.is_file() == True:
    ax2.plot(simdata[3]["T"],-(0.4/2.45)*2.65*(simdata[3]["V[8,1]"]-simdata[3]["V["+str(Nx-7)+",1]"])/simdata[3]["I"],color='purple',marker='^',label='long c-axis') 
RVTaxes(ax2)
# ax1.set_ylim(bottom=0, top=None)
# ax2.set_ylim(bottom=-10, top=None)
fig.savefig('../../Plots/Sr327/Measurements/PaperRhoC.svg')
plt.show()
