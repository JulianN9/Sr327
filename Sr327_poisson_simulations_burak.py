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

def simulate(Vin,Vout,P):
    datac = []; datad = []; datax1 = []; dataL = []
    data = [ datac,datad,datax1,dataL]
    N = 100
    # Rx = rx(T,Nx); Ry = newry(T,Ny)

    # typestr = 'c'
    typestrlist = ['c','d','x1','L']
    iolist = [inputlist(t,Nx,Ny,Vin,Vout) for t in typestrlist]
    for Tctr in range(0,N):
        T = 300.0-Tctr*(300.-2.)/(N-1)
        print(T)
        Rx = rx(T,Nx); Ry = newry(T,Ny)/(1+(P-1)/10)
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

if __name__ == "__main__":
    Vin = 5; Vout = 6
    P = 10

    for P in range(100):
        print("Pressure: "+str(P+1))
        simulate(Vin,Vout,P+1)