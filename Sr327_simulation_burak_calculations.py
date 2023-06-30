import os
import scipy
import numpy as np
import sympy as sp
from sympy import *
from sympy.matrices import Matrix
import pandas as pd
import matplotlib.pyplot as plt
from math import floor
# from multiprocessing import Pool
# from itertools import product, repeat, permutations, chain
from Sr327_simulation import rx, newry, R_ab, newR_c

meshsize = 2
Nx = 5
Ny = 3

N = Nx*Ny
A = np.zeros((N,N))
A = sp.Matrix(A)
x = sp.symbols('x')
x21 = sp.Add(sp.Mul(-2,x),-1)
x22 = sp.Add(sp.Mul(-2,x),-2)

# Diagonal
for i in range(1,Nx-1):
    A[i+Nx*(Ny-1),i+Nx*(Ny-1)]=x21
    A[i,i]=x21
    for j in range(1,Ny-1):
        A[i+Nx*j,i+Nx*j]=x22
A[Nx-1,Nx-1]=-x-1
A[0,0]=-x-1
A[Nx-1+Nx*(Ny-1),Nx-1+Nx*(Ny-1)]=-x-1
A[Nx*(Ny-1),Nx*(Ny-1)]=-x-1
for j in range(1,Ny-1):
    A[Nx-1+Nx*j,Nx-1+Nx*j]=-x-2
    A[Nx*j,Nx*j]=-x-2

# Tri-Diagonal
for i in range(0,Nx-1):
    for j in range(0,Ny):
        A[i+Nx*j+1,i+Nx*j]=x
        A[i+Nx*j,i+Nx*j+1]=x

# Off-Diagonal
for j in range(0,Ny-1):
    for i in range(0,Nx):
        A[i+Nx*j,i+Nx*(j+1)]=1
        A[i+Nx*(j+1),i+Nx*j]=1


L = 1
B = np.zeros((N-1,N-1))
B = sp.Matrix(B)
for i in range(Nx*Ny-1):
    for j in range(Nx*Ny-1):
        if (i < L)&(j<L):
            B[i,j] = A[i,j]
        elif (i >= L)&(j<L):
            B[i,j] = A[i+1,j]
        elif (i < L)&(j >= L):
            B[i,j] = A[i,j+1]
        elif (i >= L)&(j>=L):
            B[i,j] = A[i+1,j+1]
M = sp.Matrix(B)

print(M)
print(M.det())
