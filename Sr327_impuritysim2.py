import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate
import sympy
import quadpy

delta = 10**(-30)

def freal(x,t,T,debye):
    # return ((x/debye**2) * ( (1-np.cos(x*t)) * ( 1/(np.tanh(x/(2*T))) - 1 ) +(1-np.exp(1j * x * t)))).real
    return (x/debye**2) * (1-np.cos(x*t)) * ( 1/(np.tanh(x/(2*T))) ) 

def fimag(x,t,T,debye):
    return -(x/debye**2) * np.sin( x * t)

def sigma(t,T,debye,lambdap, real=True, printv = False):
    intcalcreal, interr = integrate.quad(freal,0,debye,args=(t,T,debye))
    # print("Intcalcreal " + str(t) + ": "+str(intcalcreal))
    intcalcimag = (debye*t*np.cos(debye*t)-np.sin(debye*t))/(debye**2 * t**2)
    # intcalc = intcalcreal + 1j*intcalcimag
    # sigmacalc = ( 1j * np.pi * T**2 * t)/ ( (np.sinh(np.pi*T*t))**2 ) * np.exp(-lambdap*1j*intcalc) 
    sigmacalc = -(( np.pi * T**2 * t)/ ( (np.sinh(np.pi*T*t))**2 )) * np.sin(-lambdap*intcalcimag)  * np.exp(-lambdap*intcalcreal)
    # sigmacalc = -(( np.pi * T**2 * t)/ ( (np.sinh(np.pi*T*t))**2 )) * np.sin(-lambdap*intcalcimag) 
    if printv == True:
        print("intcalcreal: "+str(intcalcreal))
        print("sigmacalc: "+str(sigmacalc))
    if real == True:
        return -sigmacalc.real
    else:
        return sigmacalc.imag

def resistivity(T):
    calcsigma, sigmaerr = integrate.quad(sigma,-0.01,0.01,args=(T))
    print('SIGMA: '+str(el*calcsigma))
    return (el*calcsigma)**(-1)

def resistivity2(T):
    intcalc, interr = integrate.quad(f,0,debye,args=(0,T))
    print(intcalc)
    return np.exp(-lambdap*intcalc)

def resistivity3(T):
    calcsigma1, sigmaerr1 = integrate.quad(sigma,-np.inf,-0.0000000001,args=(T))
    calcsigma2, sigmaerr2 = integrate.quad(sigma,0.0000000001,np.inf,args=(T))
    print('SIGMA1: '+str(el*calcsigma1))
    print('SIGMA2: '+str(el*calcsigma2))
    return (el*(calcsigma1+calcsigma2))**(-1)

def resistivity4(T,debye,lambdap,el):
    calcsigma1, sigmaerr1 = integrate.quad(sigma,0.5*10**(-10),1,args=(T,debye,lambdap))
    print('SIGMA '+str(T)+': '+str(2*el*calcsigma1))
    return (el*(2*calcsigma1))**(-1)

def resistivity5(T,debye,lambdap,el):
    calcsigma1, sigmaerr1 = integrate.quad(sigma,10**(-10),1,args=(T,debye,lambdap))
    calcsigma2, sigmaerr1 = integrate.quad(sigma,10**(-100),10^(-10),args=(T,debye,lambdap))
    calcsigma3, sigmaerr1 = integrate.quad(sigma,10**(-1000),10^(-100),args=(T,debye,lambdap))
    print('SIGMA '+str(T)+': '+str(2*el*(calcsigma1+calcsigma2+calcsigma3)))
    return (el*(2*calcsigma1+calcsigma2+calcsigma3))**(-1)

def resistivity6(T,debye,lambdap,el):
    calcsigma1, sigmaerr1 = integrate.quad(sigma,-10,10,args=(T,debye,lambdap),points=0)
    print('SIGMA '+str(T)+': '+str(el*calcsigma1))
    return (el*calcsigma1)**(-1)

def plotf(a,b,t,which=True):
    line = np.linspace(a,b,100)
    Y = []
    if which == True:
        for i in line:
            Y.append(freal(i,t,300,140))
    else:
        for i in line:
            Y.append(fimag(i,t,300,140))
    plt.scatter(line,Y)
    plt.show()
    return 0

def plotsigma(a,b):
    line = np.linspace(a,b,100)
    Y = []
    for i in line:
        Y.append(sigma(i,100,140,18.9))
    plt.scatter(line,Y)
    plt.show()
    return 0

def plotres(P):
    debye = 140
    lambdap = 18.9
    el = 1.5
    delta = 10*(-6)
    debye = debye * P
    lambdap = lambdap*P**2
    el = el
    line = np.linspace(20,300,20)
    y = []
    for temperature in line:
        y.append(resistivity4(temperature,debye,lambdap,el))

    fig, ax = plt.subplots()
    ax.scatter(line,y, c='b')
    ax.set_xlabel('Temperature (K)',fontsize=14)
    ax.set_ylabel('Resistivity (A. U.)',fontsize=14)
    # plt.savefig("../../Plots/Sr327/resistivity"+str(P)+".png")
    plt.show()
    return 0

plotres(1)

# Pressurelist = [0.95,0.995,0.9995,1.0005,1.005,1.05]
# for pressure in Pressurelist:
#     plotres(pressure)