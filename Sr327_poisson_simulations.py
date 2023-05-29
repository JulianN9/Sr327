import os
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from math import exp, ceil
from multiprocessing import Pool
from itertools import product, repeat, permutations, chain
from Sr327_simulation import rx, newry, R_ab, newR_c

Nx = 24
Ny = 5

# Set up and invert the matrix
def poissonmatrix(Nx,Ny,rx,ry,plot=False):
    N = Nx*Ny
    A = np.zeros((N,N))
    gamma = ry/rx
    # Diagonal
    for i in range(1,Nx-1):
        A[i+Nx*(Ny-1),i+Nx*(Ny-1)]=-2*gamma-1
        A[i,i]=-2*gamma-1
        for j in range(1,Ny-1):
            A[i+Nx*j,i+Nx*j]=-2*gamma-2
    A[Nx-1,Nx-1]=-gamma-1
    A[0,0]=-gamma-1
    A[Nx-1+Nx*(Ny-1),Nx-1+Nx*(Ny-1)]=-gamma-1
    A[Nx*(Ny-1),Nx*(Ny-1)]=-gamma-1
    for j in range(1,Ny-1):
        A[Nx-1+Nx*j,Nx-1+Nx*j]=-gamma-2
        A[Nx*j,Nx*j]=-gamma-2
    for i in range(0,N-1):
            A[i+1,i]=gamma
            A[i,i+1]=gamma
    for j in range(0,Ny-1):
        for i in range(0,Nx):
            A[i+Nx*j,i+Nx*(j+1)]=1
            A[i+Nx*(j+1),i+Nx*j]=1
    
    Ainv=np.linalg.inv(A)   

    if plot==True:
        fig = plt.figure(figsize=(12,4))
        plt.subplot(121)
        plt.imshow(A,interpolation='none')
        clb=plt.colorbar()
        clb.set_label('Matrix elements values')
        plt.title('Matrix A ',fontsize=24)
        plt.subplot(122)
        plt.imshow(Ainv,interpolation='none')
        clb=plt.colorbar()
        clb.set_label('Matrix elements values')
        plt.title(r'Matrix $A^{-1}$ ',fontsize=24)

        fig.tight_layout()
        plt.show()
    return Ainv

# def boundaryconditions(Nx,Ny,V_in,V_out):
#     Q = np.zeros(Nx*Ny)
#     Q[V_in] = 1
#     Q[-V_out] = 1
#     return Q

def chargecalculator(Ainv,Vin,Vout,typestr = 'c'):
    if typestr == 'c':
        M = np.matrix([[Ainv[Vin,Vin],Ainv[Vin,-(Nx-Vout)]],[Ainv[-(Nx-Vout),Vin],Ainv[-(Nx-Vout),-(Nx-Vout)]]])
    elif typestr == 'd':
        M = np.matrix([[Ainv[Vin,Vin],Ainv[Vin,-Vout]],[Ainv[-Vout,Vin],Ainv[-Vout,-Vout]]])
    elif typestr == 'x1':
        M = np.matrix([[Ainv[Vin,Vin],Ainv[Vin,Nx-Vin]],[Ainv[Nx-Vin,Vin],Ainv[Nx-Vin,Nx-Vin]]])
    Minv = np.linalg.inv(M)
    V = np.array([-0.5,0.5])
    Q = np.dot(Minv,V)
    return Q

def chargecalculatorL(Ainv,L):
    M = np.matrix([[ Ainv[i,j] for i in L ] for j in L])
    Minv = np.linalg.inv(M)
    V = []
    for i in range(0,int(len(L)/2)):
        V.append(-0.5)
    for i in range(0,int(len(L)/2)):
        V.append(0.5)
    V = np.array(V)
    Q = np.dot(Minv,V)
    return Q

def voltagematrix(Ainv,Q,Nx,Ny,Vin,Vout,typestr = 'c'):
    Qprime = np.zeros(Nx*Ny)
    if typestr == 'c':
        Qprime[Vin] = Q[0,0]
        Qprime[-(Nx-Vout)] = Q[0,1]
    if typestr == 'd':
        Qprime[Vin] = Q[0,0]
        Qprime[-Vout] = Q[0,1]
    if typestr == 'x1':
        Qprime[Vin] = Q[0,0]
        Qprime[Nx-Vin] = Q[0,1]    
    V = np.dot(Ainv,Qprime)
    V = V.reshape(Ny,Nx)
    return(V)

def voltagematrixL(Ainv,Q,Nx,Ny,L):
    Qprime = np.zeros(Nx*Ny)
    for count,value in enumerate(L):
        Qprime[value] = Q[0,count]
    V = np.dot(Ainv,Qprime)
    V = V.reshape(Ny,Nx)
    return(V)

# def savevoltagematrix(V,Nx,Ny,typestr,I_in,I_out,T,Rx,Ry):

if __name__ == "__main__":
    T = 2; Vin = 5; Vout = 6
    data = []; datac = []; datad = []; datax1 =[]; dataL = []
    N = 100
    # typestr = 'c'
    typestrlist = ['c','d','x1','L']
    for Tctr in range(0,N):
        T = 300.0-Tctr*(300.-2.)/(N-1)
        Rx = rx(T,Nx)/2; Ry = newry(T,Ny)
        Ainv = poissonmatrix(Nx,Ny,Rx,Ry)
        for typestr in typestrlist:
            if typestr == 'L':
                AinvL = poissonmatrix(Nx,Ny,Ry,Rx)
                L1 = [i*Nx for i in range(1,Ny-1)]
                L2 = [(i+1)*Nx-1 for i in range(1,Ny-1)]
                L = L1+L2
                Q = chargecalculatorL(AinvL,L)
                V = voltagematrixL(AinvL,Q,Nx,Ny,L)
                Vlist = [T,Rx,Ry,np.abs((V[1,0]-V[1,Nx-1])/Ry)]
                for i in range(1,Nx+1):
                    for j in range(1,Ny+1):
                        Vlist.append(V[j-1,i-1])
                dataL.append(Vlist)
            else:
                Q = chargecalculator(Ainv,Vin,Vout,typestr)/Ry
                # print(Q[0,0])
                V = voltagematrix(Ainv,Q,Nx,Ny,Vin,Vout,typestr)
                # savevoltagematrix(V,Nx,Ny,'c',Vin,Vout,T,Rx,Ry)
                # plt.contourf(V,cmap='RdGy')
                # plt.show()
                Vlist = [T,Rx,Ry,np.abs((V[0,Vin]-V[1,Vin])/Ry)]
                for i in range(1,Nx+1):
                    for j in range(1,Ny+1):
                        Vlist.append(V[j-1,i-1])
                if typestr == 'c':
                    datac.append(Vlist)
                elif typestr == 'd':
                    datad.append(Vlist)
                elif typestr == 'x1':
                    datax1.append(Vlist)

    headerlist = ['T','rx','ry','I'] # Header for the output file
    for i in range(1,Nx+1):
        for j in range(1,Ny+1):
            headerlist.append('V['+str(i)+','+str(j)+']')

    path = '../../Data/Sr327_ImplicitSimulator/'+str(Nx)+'x'+str(Ny)+'DAT/'
    
    for typestr in typestrlist:
        if typestr == 'c':
            df = pd.DataFrame(datac) # Making a dataframe from the data
        elif typestr == 'd':
            df = pd.DataFrame(datad) # Making a dataframe from the data
        elif typestr == 'x1':
            df = pd.DataFrame(datax1) # Making a dataframe from the data
        elif typestr == 'L':
            df = pd.DataFrame(dataL) # Making a dataframe from the data
        # Output path name
        if os.path.exists(path) == False:
            os.mkdir(path)
        df.to_csv(path+'Sr327_'+typestr+'_'+str(Vin)+'_to_'+str(Vout)+'.dat', index=False, header=headerlist)