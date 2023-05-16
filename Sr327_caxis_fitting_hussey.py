import numpy as np
import matplotlib.pyplot as plt
from math import exp, ceil, floor
from scipy import integrate
from scipy import interpolate
from scipy.optimize import curve_fit
from Sr327 import rhofactor, loadqualifierdata
import sympy
import pandas as pd

# If this is right then the effective result is two different paths that compete for transport. The law here is resistors in parallel. This indicates that the general resistivity is like rho a-c, but then there is some other method of transport.

j_path = '../../Data/Sr327/Julian_Tom_data/'
d_path = '../../Data/Sr327/Di_Tom_data/'
A = 1260*(23)**3*(10**(-3))**2
L = 0.46*10**(-3)

def R_ab(T): 
  f = 1/(np.exp((T-20)/10) + 1)
  return (1.7+0.03*T*T)*f + 0.68*T*(1.0-f)

# def R_ab(T): 
#   return T**2

def Gamma(T,gamma):
    return gamma*R_ab(T)

def ATGamma(T,gamma,C):
    return gamma*R_ab(T)+C

def AT(T,A,B):
    return A*np.log(B/T)

def ATexp(T,A,B,C):
    return A*np.exp(C*T)+B

def ATlin(T,A,B):
    return A*T+B

def fit(T,gamma,A,B):
    rho = 1/(1/Gamma(T,gamma)+1/ATGamma(T,A,B))
    return rho

#Reading Data
# path = '../../Data/Julian_Sr327_Tom_data_2022/isolation_c_cooling_down_2/Data.csv'
# header = ['Index','Temperature','SignalX','SignalY','Capacitancex','Capacitancey']
# data = pd.read_csv(path, skiprows=1, names=header)
# rhoc = np.sqrt(data[header[2]]**2+data[header[3]])
# rhoc = rhoc.to_numpy()
# T = data[header[1]]

data = loadqualifierdata()
rhoc_raw = data[4][0]['SignalR']*data[4][2]
T_raw = data[4][0]['Temperature']
T = np.linspace(4.5,291,100)
R_i = interpolate.interp1d(T_raw,rhoc_raw)
rhoc = R_i(T)

guess = [1/1.5,0.01,9]
params, cov = curve_fit(fit,T,rhoc, p0=guess)
fitresults = fit(T,params[0],params[1],params[2])
# fitresults = fitresults.to_numpy()

print(params)
print(cov)

# T = np.linspace(1,300)

fig = plt.figure(figsize=(10,8))
a0, a1 = fig.subplots(2,1, sharex='col', gridspec_kw={'height_ratios': [3, 1]})

a0.plot(T,rhoc,ls=':')
a0.plot(T,fit(T,params[0],params[1],params[2]))
a0.plot(T,Gamma(T,params[0]),c='k',ls='--')
a0.plot(T,ATGamma(T,params[1],params[2]),c='k',ls='--')
a0.plot(T,Gamma(T,guess[0]),c='k',ls='-.')
a0.plot(T,ATGamma(T,guess[1],guess[2]),c='k',ls='-.')
a0.set_ylim(0,13)
a0.set_xlim(0,300)
a0.tick_params(axis='x', labelsize=16)
a0.tick_params(axis='y', labelsize=16)
a0.set_ylabel('Resistivity (m$\Omega$cm)',fontsize=24)
# plt.show()

# Residuals:
residual = rhoc - fitresults
a1.plot(T,residual)
a1.axhline(c='k')
a1.set_xlim(0,300)
a1.tick_params(axis='x', labelsize=16)
a1.tick_params(axis='y', labelsize=16)
a1.set_xlabel('Temperature (K)',fontsize=24)
fig.savefig('../../Plots/Sr327/HusseyFitting.svg')
plt.show()