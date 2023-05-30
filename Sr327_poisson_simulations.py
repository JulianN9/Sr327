import os
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from math import exp, ceil, floor
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

def chargecalculator(Ainv,L):
    M = np.matrix([[ Ainv[j,i] for i in L ] for j in L])
    Minv = np.linalg.inv(M)
    V = []
    for i in range(0,int(len(L)/2)):
        V.append(-0.5)
    for i in range(0,int(len(L)/2)):
        V.append(0.5)
    V = np.array(V)
    Q = np.dot(Minv,V)
    return Q

def voltagematrix(Ainv,Q,Nx,Ny,L):
    Qprime = np.zeros(Nx*Ny)
    for count,value in enumerate(L):
        Qprime[value] = Q[0,count]
    V = np.dot(Ainv,Qprime)
    V = V.reshape(Ny,Nx)
    return(V)

def inputlist(typestr,Nx,Ny,Vin,Vout):
    if typestr == 'c':
        L = [Vin,-(Nx-Vout)]
    elif typestr == 'd':
        L = [Vin,-Vout]
    elif typestr == 'x1':
        L = [Vin,Nx-Vin]
    elif typestr == 'L':
        L1 = [i*Nx for i in range(1,Ny-1)]
        L2 = [(i+1)*Nx-1 for i in range(1,Ny-1)]
        L = L1+L2
    else:
        print('ERROR')
    return L

# def savevoltagematrix(V,Nx,Ny,typestr,I_in,I_out,T,Rx,Ry):

if __name__ == "__main__":
    T = 2; Vin = 5; Vout = 6
    datac = []; datad = []; datax1 = []; dataL = []
    data = [ datac,datad,datax1,dataL]
    N = 100
    # typestr = 'c'
    typestrlist = ['c','d','x1','L']
    iolist = [inputlist(t,Nx,Ny,Vin,Vout) for t in typestrlist]
    for Tctr in range(0,N):
        T = 300.0-Tctr*(300.-2.)/(N-1)
        Rx = rx(T,Nx)/2; Ry = newry(T,Ny)
        Ainv = poissonmatrix(Nx,Ny,Rx,Ry)
        AinvL = poissonmatrix(Nx,Ny,Ry,Rx)
        Alist = [Ainv, Ainv, Ainv, AinvL]
        for count, typestr in enumerate(typestrlist):
            Q = chargecalculator(Alist[count],iolist[count])
            V = voltagematrix(Alist[count],Q,Nx,Ny,iolist[count])

            # savevoltagematrix(V,Nx,Ny,'c',Vin,Vout,T,Rx,Ry)
            # plt.contourf(V,cmap='RdGy')
            # plt.show()

            if (typestr == 'c')|(typestr == 'd'): 
                Vlist = [T,Rx,Ry,np.abs((V[0,iolist[count][0]]-V[1,iolist[count][0]])/Ry)]
            elif typestr == 'x1': 
                Vlist = [T,Rx,Ry,np.abs((V[0,iolist[count][0]]-V[0,iolist[count][0]+1])/Ry)]
            elif typestr == 'L':
                Vlist = [T,Rx,Ry,np.abs((V[floor(Ny/2),0]-V[floor(Ny/2),1])/Ry)]
            for i in range(1,Nx+1):
                for j in range(1,Ny+1):
                    Vlist.append(V[j-1,i-1])
            data[count].append(Vlist)

    headerlist = ['T','rx','ry','I'] # Header for the output file
    for i in range(1,Nx+1):
        for j in range(1,Ny+1):
            headerlist.append('V['+str(i)+','+str(j)+']')

    path = '../../Data/Sr327_ImplicitSimulator/'+str(Nx)+'x'+str(Ny)+'DAT/'
    
    for count, typestr in enumerate(typestrlist):
        df = pd.DataFrame(data[count]) # Making a dataframe from the data
        # Output path name
        if os.path.exists(path) == False:
            os.mkdir(path)
        df.to_csv(path+'Sr327_'+typestr+'_'+str(Vin)+'_to_'+str(Vout)+'.dat', index=False, header=headerlist)