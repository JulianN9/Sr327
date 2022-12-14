import argparse
from itertools import product, repeat
from multiprocessing import Pool
from timeit import repeat
from math import exp, ceil
import numpy as np

def R_c(T):  # in mu-Ohm-cm (hence factor of 1000.0)
  f = 1.0/(exp((T-10.0)/10.) + 1.)
  g = 1.0/(exp((T-40.0)/20.) + 1.)
  return 1000.0*(f*(1.0+0.02*T*T) + 0.08*T*g*(1.-f) + (8.0+0.0005*T)*(1.-g)*(1.-f))

def R_ab(T): 
  f = 1.0/(exp((T-20.)/10.) + 1.0)
  return (1.7+0.03*T*T)*f + 0.68*T*(1.0-f)

Nx = 12
Ny = 3
I_in = 3
I_out = 4
T = 300


V = np.matrix(np.zeros([Nx+2,Ny+2]))
dQ = np.matrix(np.zeros([Nx+2,Ny+2]))
V[I_in+1,1] = -0.5 
V[I_out+1,-2] = 0.5  # 1
rx = 200*R_ab(T)/(Nx-1)
ry = 20*R_c(T)/(Ny-1)
V[1:Nx+1,1:Ny+1] += 10*( (-2*V[1:Nx+1,1:Ny+1] + V[0:Nx,1:Ny+1] + V[2:Nx+2,1:Ny+1])/rx + (-2*V[1:Nx+1,1:Ny+1] + V[1:Nx+1,0:Ny] + V[1:Nx+1,2:Ny+2])/ry  ) 
print(dQ)
print(V)