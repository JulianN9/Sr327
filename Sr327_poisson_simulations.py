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

# Set up and invert the matrix
def poissonmatrix(Nx,Ny,rx,ry,L,plot=False):
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

    # Tri-Diagonal
    for i in range(0,Nx-1):
        for j in range(0,Ny):
            A[i+Nx*j+1,i+Nx*j]=gamma
            A[i+Nx*j,i+Nx*j+1]=gamma
    
    # Off-Diagonal
    for j in range(0,Ny-1):
        for i in range(0,Nx):
            A[i+Nx*j,i+Nx*(j+1)]=1
            A[i+Nx*(j+1),i+Nx*j]=1

    for pin in L:
        A[pin,:] = 0
        A[pin,pin] = 1

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
    return A, Ainv

# Set up and decompose the matrix
def poissonmatrixLU(Nx,Ny,rx,ry,plot=False):
    N = Nx*Ny
    A = np.zeros((N,N),dtype=np.longdouble)
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

    # Tri-Diagonal
    for i in range(0,Nx-1):
        for j in range(0,Ny):
            A[i+Nx*j+1,i+Nx*j]=gamma
            A[i+Nx*j,i+Nx*j+1]=gamma
    
    # Off-Diagonal
    for j in range(0,Ny-1):
        for i in range(0,Nx):
            A[i+Nx*j,i+Nx*(j+1)]=1
            A[i+Nx*(j+1),i+Nx*j]=1
    
    P, L, U = scipy.linalg.lu(A)
    L = scipy.linalg.inv(L)
    U = scipy.linalg.inv(U)

    if plot==True:
        fig = plt.figure(figsize=(12,4))
        plt.subplot(131)
        plt.imshow(A,interpolation='none')
        clb=plt.colorbar()
        clb.set_label('Matrix elements values')
        plt.title('Matrix A ',fontsize=24)
        plt.subplot(132)
        plt.imshow(U,interpolation='none')
        clb=plt.colorbar()
        clb.set_label('Matrix elements values')
        plt.title(r'Matrix U ',fontsize=24)
        plt.subplot(133)
        plt.imshow(L,interpolation='none')
        clb=plt.colorbar()
        clb.set_label('Matrix elements values')
        plt.title(r'Matrix L ',fontsize=24)

        fig.tight_layout()
        plt.show()
    return A,L,U

def poissonmatrixSCHUR(Nx,Ny,rx,ry,plot=False):
    N = Nx*Ny
    A = np.zeros((N,N),dtype=np.longdouble)
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

    # Tri-Diagonal
    for i in range(0,Nx-1):
        for j in range(0,Ny):
            A[i+Nx*j+1,i+Nx*j]=gamma
            A[i+Nx*j,i+Nx*j+1]=gamma
    
    # Off-Diagonal
    for j in range(0,Ny-1):
        for i in range(0,Nx):
            A[i+Nx*j,i+Nx*(j+1)]=1
            A[i+Nx*(j+1),i+Nx*j]=1
    
    T, Z = scipy.linalg.schur(A)
    T = scipy.linalg.inv(T)

    if plot==True:
        fig = plt.figure(figsize=(12,4))
        plt.subplot(131)
        plt.imshow(A,interpolation='none')
        clb=plt.colorbar()
        clb.set_label('Matrix elements values')
        plt.title('Matrix A ',fontsize=24)
        plt.subplot(132)
        plt.imshow(T,interpolation='none')
        clb=plt.colorbar()
        clb.set_label('Matrix elements values')
        plt.title(r'Matrix T ',fontsize=24)
        plt.subplot(133)
        plt.imshow(Z,interpolation='none')
        clb=plt.colorbar()
        clb.set_label('Matrix elements values')
        plt.title(r'Matrix Z ',fontsize=24)

        fig.tight_layout()
        plt.show()
    return A,T,Z

# Set up and invert the matrix
def poissonmatrixPSEUDO(Nx,Ny,rx,ry,plot=False):
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

    # Tri-Diagonal
    for i in range(0,Nx-1):
        for j in range(0,Ny):
            A[i+Nx*j+1,i+Nx*j]=gamma
            A[i+Nx*j,i+Nx*j+1]=gamma
    
    # Off-Diagonal
    for j in range(0,Ny-1):
        for i in range(0,Nx):
            A[i+Nx*j,i+Nx*(j+1)]=1
            A[i+Nx*(j+1),i+Nx*j]=1
    
    Ainv=np.linalg.pinv(A)   

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
    return A, Ainv

def chargecalculator(Linv,Uinv,L,N,Rc):
    Nq = len(L)
    Selector = np.zeros((N,Nq))
    for count,value in enumerate(L):
        Selector[value,count] = 1
    SelectorT = np.transpose(Selector)
    M = np.matmul(SelectorT,np.matmul(Uinv,np.matmul(Linv,Selector)))
    V = []
    for i in range(0,int(Nq/2)):
        V.append(-0.5*scale/Rc)
    for i in range(0,int(Nq/2)):
        V.append(0.5*scale/Rc)
    V = np.array(V)
    lu, piv = scipy.linalg.lu_factor(M)
    Q = scipy.linalg.lu_solve((lu, piv),V)
    return Q

def chargecalculatorPSEUDO(Ainv,L,N,Rc):
    Nq = len(L)
    Selector = np.zeros((N,Nq))
    for count,value in enumerate(L):
        Selector[value,count] = 1
    SelectorT = np.transpose(Selector)
    M = np.matmul(SelectorT,np.matmul(Ainv,Selector))
    V = []
    for i in range(0,int(Nq/2)):
        V.append(-0.5*scale/Rc)
    for i in range(0,int(Nq/2)):
        V.append(0.5*scale/Rc)
    V = np.array(V)
    lu, piv = scipy.linalg.lu_factor(M)
    Q = scipy.linalg.lu_solve((lu, piv),V)
    return Q

def chargecalculatorTZ(Tinv,Z,L,N,Rc):
    Zinv = np.transpose(Z)
    Nq = len(L)
    Selector = np.zeros((N,Nq))
    for count,value in enumerate(L):
        Selector[value,count] = 1
    SelectorT = np.transpose(Selector)
    M = np.matmul(SelectorT,np.matmul(Z,np.matmul(Tinv,np.matmul(Zinv,Selector))))
    V = []
    for i in range(0,int(Nq/2)):
        V.append(-0.5*scale/Rc)
    for i in range(0,int(Nq/2)):
        V.append(0.5*scale/Rc)
    V = np.array(V)
    lu, piv = scipy.linalg.lu_factor(M)
    Q = scipy.linalg.lu_solve((lu, piv),V)
    return Q

def voltagematrix(Linv,Uinv,L,Q,Nx,Ny):
    Qprime = np.zeros(Nx*Ny)
    for count,value in enumerate(L):
        Qprime[value] = Q[count]
    b = np.dot(Linv,Qprime)
    V = np.dot(Uinv,b)
    # print(type(V))
    # V = V.reshape(Ny,Nx)
    return(V)

def voltagematrixPSEUDO(Ainv,L,Q,Nx,Ny):
    Qprime = np.zeros(Nx*Ny)
    for count,value in enumerate(L):
        Qprime[value] = Q[count]
    V = np.matmul(Ainv,Qprime)
    # print(type(V))
    # V = V.reshape(Ny,Nx)
    return(V)


def voltagematrixTZ(Tinv,Z,L,Q,Nx,Ny):
    Zinv = np.transpose(Z)
    Qprime = np.zeros(Nx*Ny)
    for count,value in enumerate(L):
        Qprime[value] = Q[count]
    b = np.matmul(Zinv,Qprime)
    c = np.matmul(Tinv,b)
    V = np.matmul(Z,c)
    # print(type(V))
    # V = V.reshape(Ny,Nx)
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

def convergence_speed(T):
   cs = 200/(1+np.exp(-(1/40)*(T-150)))
   return cs

# Convergence functions which give the final values of voltage across the lattice, c is for c-axis setup, d is for diagonal, x1 and x2 are for in-plane, long and short
def converge(V, A, L, T, RX, RY):
    # cs = convergence_speed(T) # Get convergence speed
    err = 1.0; ctr = 0 # This defines the change between the guess pin and it's neighbors and the count, these are used to determine convergence.
    dQerr = 0; dQperr = 0
    while ((err > 10**(-7 - floor(len(L)/2)))|(ctr<100))&(ctr<200000):
        # csp = (cs/(1+np.exp(-(1/500)*(10000-ctr))))+cs
        for count, value in enumerate(L):
            if count < floor(len(L)/2):
                V[value] = -0.5
            else:
                V[value] = 0.5
        dQ = np.matmul(A,V)
        V = dQ + V

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

if __name__ == "__main__":
    T = 2; Vin = 5; Vout = 6
    P = 1
    datac = []; datad = []; datax1 = []; dataL = []
    data = [ datac,datad,datax1,dataL]
    N = 100
    # Rx = rx(T,Nx)/2; Ry = newry(T,Ny)
    # Ainv = poissonmatrix(Nx,Ny,1,2,True)

    # typestr = 'c'
    typestrlist = ['c','d','x1','L']
    iolist = [inputlist(t,Nx,Ny,Vin,Vout) for t in typestrlist]
    for Tctr in range(0,N):
        T = 300.0-Tctr*(300.-2.)/(N-1)
        print(T)
        Rx = rx(T,Nx)/2; Ry = newry(T,Ny)/P
        Rxlist = [Rx, Rx, Rx, Ry]
        Rylist = [Ry, Ry, Ry, Rx]
        # Alist = [A, A, A, AL]
        # Ainvlist = [Ainv, Ainv, Ainv, AinvL]
        Rc = [Ry,Ry,Ry,Rx]
        for count, typestr in enumerate(typestrlist):
            A, Ainv = poissonmatrix(Nx,Ny,Rxlist[count],Rylist[count],iolist[count],False)
            Q = chargecalculatorPSEUDO(Ainv,iolist[count],Nx*Ny,Rc[count])
            V = voltagematrixPSEUDO(Ainv,iolist[count],Q,Nx,Ny)
            dQ = np.matmul(A,V)
            V = V.reshape(Ny,Nx)
            # plt.contourf(V,cmap='RdGy')
            # plt.show()

            npins = floor(len(iolist[count])/2)
            Vlist = [T,Rx,Ry,np.sum([dQ[iolist[count][i]] for i in range(npins)])/npins]
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