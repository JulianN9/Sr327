import os
import scipy
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from math import floor
# from multiprocessing import Pool
# from itertools import product, repeat, permutations, chain
from Sr327_simulation import rx, newry, R_ab, newR_c

meshsize = 2
Nx = 12*meshsize
Ny = 2*meshsize+1
scale = 10**(0)

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

def convergence_speed(T):
    cs = 200/(1+np.exp(-(1/40)*(T-150)))
    # cs = 100
    return cs

# Convergence functions which give the final values of voltage across the lattice, c is for c-axis setup, d is for diagonal, x1 and x2 are for in-plane, long and short
def converge(V, A, L, T, RX, RY):
    # if floor(len(L)/2) > 1:
    #     cs = 10 # Get convergence speed
    # else:
    cs = convergence_speed(T) # Get convergence speed
    # rxt = rx(T,Nx); ryt = newry(T,Ny) # Get strength of resistors
    err = 1.0; ctr = 0 # This defines the change between the guess pin and it's neighbors and the count, these are used to determine convergence.
    dQerr = 0; dQperr = 0
    while ((err > 10**(-11 - floor(len(L)/2)))|(ctr<1000))&(ctr<200000):
        csp = (cs/(1+np.exp(-(1/500)*(10000-ctr))))+cs
        for count, value in enumerate(L):
            if count < floor(len(L)/2):
                V[value] = -0.5
            else:
                V[value] = 0.5
        dQ = np.matmul(A,V)
        V = csp*dQ + V

        dQperr = dQerr
        dQerr = np.sum(np.abs(dQ))/(Nx*Ny)
        err = dQperr - dQerr
        # err = dQerr
        # print(err)

        ctr += 1
    Vlist = [T,RX,RY,np.sum([dQ[L[i]] for i in range(floor(len(L)/2))])/floor(len(L)/2)] # Setting the output lists
    V = V.reshape(Ny,Nx)
    for i in range(1,Nx+1):
        for j in range(1,Ny+1):
            Vlist.append(V[j-1,i-1])
    return Vlist

# def savevoltagematrix(V,Nx,Ny,typestr,I_in,I_out,T,Rx,Ry):

def simulate(Vin, Vout, P):
    datac = []; datad = []; datax1 = []; dataL = []
    data = [ datac,datad,datax1,dataL]
    N = 100

    # typestr = 'c'
    typestrlist = ['c','d','x1','L']
    # typestrlist = ['c','d','x1']
    iolist = [inputlist(t,Nx,Ny,Vin,Vout) for t in typestrlist]
    for Tctr in range(0,N):
        T = 300.0-Tctr*(300.-2.)/(N-1)
        print(T)
        Rx = rx(T,Nx); Ry = newry(T,Ny)/(1+(P-1)/10)
        A = poissonmatrix(Nx,Ny,Rx,Ry)
        AL = poissonmatrix(Nx,Ny,Ry,Rx)
        Alist = [A, A, A, AL]
        Rc = [Ry,Ry,Ry,Rx]
        for count, typestr in enumerate(typestrlist):
            V = np.zeros(Nx*Ny)
            data[count].append(converge(V, Alist[count], iolist[count], T, Rx, Ry))


    headerlist = ['T','rx','ry','I'] # Header for the output file
    for i in range(1,Nx+1):
        for j in range(1,Ny+1):
            headerlist.append('V['+str(i)+','+str(j)+']')

    if P!=0:
        path = '../../Data/Sr327_ImplicitSimulator/'+str(P)+'P'+str(Nx)+'x'+str(Ny)+'DAT/'
    else:
        path = '../../Data/Sr327_ImplicitSimulator/'+str(Nx)+'x'+str(Ny)+'DAT/'

    for count, typestr in enumerate(typestrlist):
        df = pd.DataFrame(data[count]) # Making a dataframe from the data
        # Output path name
        if os.path.exists(path) == False:
            os.mkdir(path)
        df.to_csv(path+'Sr327_'+typestr+'_'+str(Vin)+'_to_'+str(Vout)+'.dat', index=False, header=headerlist)
    return 0

if __name__ == "__main__":
    Vin = 5; Vout = 6 
    for P in range(10):
        simulate(Vin,Vout,P+1)
