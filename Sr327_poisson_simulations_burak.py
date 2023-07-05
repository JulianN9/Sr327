import os
import scipy
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from math import floor, exp
# from multiprocessing import Pool
# from itertools import product, repeat, permutations, chain

# Defining Resitivities:
# c-axis resistivity
def pfactor(p):
    return (1+(p-1)/10)

def R_c(T,p=1):  # in mu-Ohm-cm (hence factor of 1000.0)
    f = 1.0/(exp((T-10.0)/10.) + 1.)
    g = 1.0/(exp((T-40.0)/20.) + 1.)
    pf = pfactor(p)
    return 1000.0*(f*(1.0+0.02*T*T) + 0.08*T*g*(1.-f) + (8.0+0.0005*T)*(1.-g)*(1.-f))

def newR_c(T,p=1,alpha=1,beta=1,gamma=1):
    # [A,B,C,D,E,F,G,H] = [3.04684366e-17,10.6253646,37.8301757,12.6081855,8.31419655e-02,7.26861931e-17,7.49548715,1.51364092e-02]
    [A,B,C,D,E,F,G,H] = [ 10, 10, 40, 20,  0.022, 0.1, 7.5, 0.0005*30]
    f = 1.0/(np.exp((T-A)/B) + 1.)
    g = 1.0/(np.exp((T-C)/D) + 1.)
    pf = pfactor(p)
    return 1000.0*(f*(1.0+E*T*T)/pow(p,alpha) + F*T*g*(1.-f)/pow(p,beta) + (G+H*T)*(1.-g)*(1.-f)/pow(p,gamma))

# a-b plane resistivity
def R_ab(T,p=1): 
    f = 1.0/(exp((T-20.)/10.) + 1.0)
    pf = pfactor(p)
    return (1.7+0.03*T*T)*f + 0.68*T*(1.0-f)

# Dividing by size of lattice to get strength of resistors, 160 and 20 are arbitrary scaling factors
def rx(T,Nx,p=1):
    return 300*R_ab(T,p)/(Nx-1) #For summer 2023 data
    # return R_ab(T)/(Nx-1)

def ry(T,Ny,p=1):
    return 20*R_c(T,p)/(Ny-1) #For summer 2023 data
    # return R_c(T)/(Ny-1)

def newry(T,Ny,p=1,alpha=1,beta=1,gamma=1):
    return 20*newR_c(T,p,alpha,beta,gamma)/(Ny-1)

# Set up the matrix
def poissonmatrix(Nx,Ny,rx,ry,plot=False):
    N = Nx*Ny
    A = np.zeros((N,N))
    gamma = 1/rx
    delta = 1/ry

    # Diagonal
    for i in range(1,Nx-1):
        A[i+Nx*(Ny-1),i+Nx*(Ny-1)]=-2*gamma-delta
        A[i,i]=-2*gamma-delta
        for j in range(1,Ny-1):
            A[i+Nx*j,i+Nx*j]=-2*gamma-2*delta
    A[Nx-1,Nx-1]=-gamma-delta
    A[0,0]=-gamma-delta
    A[Nx-1+Nx*(Ny-1),Nx-1+Nx*(Ny-1)]=-gamma-delta
    A[Nx*(Ny-1),Nx*(Ny-1)]=-gamma-delta
    for j in range(1,Ny-1):
        A[Nx-1+Nx*j,Nx-1+Nx*j]=-gamma-2*delta
        A[Nx*j,Nx*j]=-gamma-2*delta

    # Tri-Diagonal
    for i in range(0,Nx-1):
        for j in range(0,Ny):
            A[i+Nx*j+1,i+Nx*j]=gamma
            A[i+Nx*j,i+Nx*j+1]=gamma
    
    # Off-Diagonal
    for j in range(0,Ny-1):
        for i in range(0,Nx):
            A[i+Nx*j,i+Nx*(j+1)]=delta
            A[i+Nx*(j+1),i+Nx*j]=delta
    
    if plot==True:
        fig = plt.figure(figsize=(12,4))
        plt.subplot(111)
        plt.imshow(A,interpolation='none')
        clb=plt.colorbar()
        clb.set_label('Matrix elements values')
        plt.title('Matrix A ',fontsize=24)

        fig.tight_layout()
        plt.show()
    return A

# Calculating the Voltage
def voltagematrix(A,L,Nx,Ny):
    V = np.zeros(Nx*Ny)
    for count,value in enumerate(L):
        if count < floor(len(L)/2):
            V[value] = -0.5
        else:
            V[value] = 0.5

    b = np.matmul(A,V)
    Ar = A
    b = np.delete(b,L,0)
    Ar = np.delete(np.delete(Ar,L,0),L,1)

    Arinv = np.linalg.inv(Ar)
    Vr = np.matmul(Arinv,b)
    rcount = 0
    for i in range(Nx*Ny):
        if (not(i in L))&(not(-Nx*Ny+i in L)):
            V[i] = -Vr[i-rcount]
        else:
            rcount = rcount + 1
    return(V)

# Defining the input/output pins
def inputlist(typestr,Nx,Ny,Vin,Vout):
    if typestr == 'c':
        L = [Vin,-(Nx-Vout)]
    elif typestr == 'd':
        L = [Vin,-Vout-1]
    elif typestr == 'x1':
        L = [Vin,Nx-1-Vin]
    elif typestr == 'L':
        L1 = [i*Nx for i in range(1,Ny-1)]
        L2 = [(i+1)*Nx-1 for i in range(1,Ny-1)]
        L = L1+L2
    else:
        print('ERROR')
    return L

# Simulated for a given set of input/output pins and pressure
def simulate(Vin,Vout,P,alpha=1,beta=1,gamma=1,save=True,verbose = False,name = ''):
    meshsize = 2
    Nx = 12*meshsize
    Ny = 2*meshsize+1
    datac = []; datad = []; datax1 = []; dataL = []
    data = [ datac,datad,datax1,dataL]
    N = 100
    # Rx = rx(T,Nx); Ry = newry(T,Ny)

    # typestr = 'c'
    typestrlist = ['c','d','x1','L']
    iolist = [inputlist(t,Nx,Ny,Vin,Vout) for t in typestrlist]
    for Tctr in range(0,N):
        T = 300.0-Tctr*(300.-2.)/(N-1)
        if verbose == True:
            print("T: "+str(T))
        Rx = rx(T,Nx); Ry = newry(T,Ny,P,alpha,beta,gamma)
        A = poissonmatrix(Nx,Ny,Rx,Ry,False)
        AL = poissonmatrix(Nx,Ny,Ry,Rx)
        Alist = [A, A, A, AL]
        Rc = [Ry,Ry,Ry,Rx]
        for count, typestr in enumerate(typestrlist):
            V = voltagematrix(Alist[count],iolist[count],Nx,Ny)

            # dQ = [ np.matmul(A,V)[i] for i in iolist[count] ]
            dQ = np.matmul(A,V)

            V = V.reshape(Ny,Nx)

            if typestr == 'L':
                Vlist = [T,Rx,Ry,np.abs((V[floor(Ny/2),0]-V[floor(Ny/2),1])/Ry)]
            else:
                Vlist = [T,Rx,Ry,np.sum(np.abs(dQ[iolist[count]]))/len(iolist[count])]
            for i in range(1,Nx+1):
                for j in range(1,Ny+1):
                    Vlist.append(V[j-1,i-1])
            data[count].append(Vlist)
            # data[count].append(converge(V,(1/Ry)*Alist[count],iolist[count],T,Rx,Ry))

    headerlist = ['T','rx','ry','I'] # Header for the output file
    for i in range(1,Nx+1):
        for j in range(1,Ny+1):
            headerlist.append('V['+str(i)+','+str(j)+']')

    if name =='':
        if P!=0:
            name = str(int(P))+'P'+str(Nx)+'x'+str(Ny)+'DAT/'
        else:
            name = str(Nx)+'x'+str(Ny)+'DAT/'
    path = '../../Data/Sr327_ImplicitSimulator/'+name

    dflist = []
    for count, typestr in enumerate(typestrlist):
        df = pd.DataFrame(data[count],columns=headerlist) # Making a dataframe from the data
        # Output path name
        if save == True:
            if os.path.exists(path) == False:
                os.mkdir(path)
            df.to_csv(path+'Sr327_'+typestr+'_'+str(Vin)+'_to_'+str(Vout)+'.dat', index=False, header=headerlist)
        dflist.append([df,typestr])
    return dflist

if __name__ == "__main__":
    Vin = 5; Vout = 6
    P = 10

    for P in range(100):
        print("Pressure: "+str(P+1))
        simulate(Vin,Vout,P+1,save=True,verbose=True)