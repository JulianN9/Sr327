import numpy as np
import matplotlib.pyplot as plt
from math import exp, ceil, floor
from scipy import integrate
from scipy.optimize import curve_fit
import sympy
import pandas as pd

# If this is right then the effective result is two different paths that compete for transport. The law here is resistors in parallel. This indicates that the general resistivity is like rho a-c, but then there is some other method of transport.

def R_ab(T): 
  f = 1/(np.exp((T-20)/10) + 1)
  return (1.7+0.03*T*T)*f + 0.68*T*(1.0-f)

# def R_ab(T): 
#   return T**2

def Gamma(T,gamma):
    return gamma*R_ab(T)

def AT(T,A,B):
    return A*np.log(B/T)

def ATexp(T,A,B,C):
    return A*np.exp(C*T)+B

def ATlin(T,A,B):
    return A*T+B

def fit(T,gamma,A,B):
    rho = 1/(1/Gamma(T,gamma)+1/AT(T,A,B))
    return rho

#Reading Data
path = '../../Data/Julian_Sr327_Tom_data_2022/isolation_c_cooling_down_2/Data.csv'
header = ['Index','Temperature','SignalX','SignalY','Capacitancex','Capacitancey']
data = pd.read_csv(path, skiprows=1, names=header)
rhoc = np.sqrt(data[header[2]]**2+data[header[3]])
rhoc = rhoc.to_numpy()
T = data[header[1]]

guess = [0.5,2000,1000]
params, cov = curve_fit(fit,T,rhoc, p0=guess)
fitresults = fit(T,params[0],params[1],params[2])
fitresults = fitresults.to_numpy()

print(params)
print(cov)

# T = np.linspace(1,300)

fig, (a0, a1) = plt.subplots(2,1, sharex='col', gridspec_kw={'height_ratios': [3, 1]})

a0.plot(T,rhoc,ls=':')
a0.plot(T,fit(T,params[0],params[1],params[2]))
a0.plot(T,Gamma(T,params[0]),c='k',ls='--')
a0.plot(T,AT(T,params[1],params[2]),c='k',ls='--')
a0.plot(T,Gamma(T,guess[0]),c='k',ls='-.')
a0.plot(T,AT(T,guess[1],guess[2]),c='k',ls='-.')
a0.set_ylim(0,5)
a0.set_xlim(0,330)
# plt.show()

# Residuals:
residual = rhoc - fitresults
a1.plot(T,residual)
a1.axhline(c='k')
plt.show()