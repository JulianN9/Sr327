import argparse
import os
import numpy as np
import pandas as pd
from math import exp, ceil
from multiprocessing import Pool
from itertools import product, repeat, permutations

# c-axis resistivity
def R_c(T):  # in mu-Ohm-cm (hence factor of 1000.0)
  f = 1.0/(exp((T-10.0)/10.) + 1.)
  g = 1.0/(exp((T-40.0)/20.) + 1.)
  return 1000.0*(f*(1.0+0.02*T*T) + 0.08*T*g*(1.-f) + (8.0+0.0005*T)*(1.-g)*(1.-f))

def newR_c(T):
   # [A,B,C,D,E,F,G,H] = [3.04684366e-17,10.6253646,37.8301757,12.6081855,8.31419655e-02,7.26861931e-17,7.49548715,1.51364092e-02]
   [A,B,C,D,E,F,G,H] = [ 10, 10, 40, 20,  0.022, 0.1, 7.5, 0.0005*30]
   f = 1.0/(np.exp((T-A)/B) + 1.)
   g = 1.0/(np.exp((T-C)/D) + 1.)
   return 1000.0*(f*(1.0+E*T*T) + F*T*g*(1.-f) + (G+H*T)*(1.-g)*(1.-f))

# a-b plane resistivity
def R_ab(T): 
  f = 1.0/(exp((T-20.)/10.) + 1.0)
  return (1.7+0.03*T*T)*f + 0.68*T*(1.0-f)

# Dividing by size of lattice to get strength of resistors, 160 and 20 are arbitrary scaling factors
def rx(T,Nx):
   return 300*R_ab(T)/(Nx-1) #For summer 2023 data
   # return R_ab(T)/(Nx-1)

def ry(T,Ny):
   return 20*R_c(T)/(Ny-1) #For summer 2023 data
   # return R_c(T)/(Ny-1)

def newry(T,Ny):
   return 20*newR_c(T)/(Ny-1)

# Convergence speed based on temperature
def convergence_speed(T):
   cs = 10
   if T < 10:
      cs = 0.5
   elif T < 50:
      cs = 1
   elif T > 250:
      cs = 20
   return cs

def convergence_speed_L(T):
   cs = 10/(1+np.exp(-(1/40)*(T-150)))
   return cs

Nx = 24 # x nodes
Ny = 5 # y nodes

# Convergence functions which give the final values of voltage across the lattice, c is for c-axis setup, d is for diagonal, x1 and x2 are for in-plane, long and short
def convergec(I_in, I_out, T, V):
   cs = convergence_speed_L(T) # Get convergence speed
   rxt = rx(T,Nx); ryt = newry(T,Ny) # Get strength of resistors
   err = 1.0; ctr = 0 # This defines the change between the guess pin and it's neighbors and the count, these are used to determine convergence.
   while ((err > 1.0e-10)|(ctr<10000))&(ctr<200000):
      csp = (cs/(1+np.exp(-(1/500)*(5000-ctr))))+cs
      V[I_in+1,1] = -0.5 # Sets input and output pins to be constant
      V[I_out+1,-2] = 0.5  # This used to be 1, and the other one used to be 1
      V[0,1:Ny+1] = V[1,1:Ny+1] # These four lines set the edge "ghost" pins to match their neighbors
      V[-1,1:Ny+1] = V[-2,1:Ny+1]
      V[1:Nx+1,0] = V[1:Nx+1,1]
      V[1:Nx+1,-1] = V[1:Nx+1,-2]
      
      # Calculating the change based on Ohm's law
      dQ = csp*( (-2*V[1:Nx+1,1:Ny+1] + V[0:Nx,1:Ny+1] + V[2:Nx+2,1:Ny+1])/rxt \
         + (-2*V[1:Nx+1,1:Ny+1] + V[1:Nx+1,0:Ny] + V[1:Nx+1,2:Ny+2])/ryt  ) 
      
      # Updating the matrix and convergence terms
      V[1:Nx+1,1:Ny+1] += dQ
      err = np.abs(np.sum(dQ[1:Nx+1,1:Ny+1])/(Nx*Ny))
      # print(err)
      ctr += 1
   Vlist = [T,rxt,ryt,dQ[I_in,0]/csp] # Setting the output lists
   for i in range(1,Nx+1):
      for j in range(1,Ny+1):
         Vlist.append(V[i,j])
   return Vlist

def converged(I_in, I_out, T, V):
   cs = convergence_speed_L(T) # Get convergence speed
   rxt = rx(T,Nx); ryt = newry(T,Ny) # Get strength of resistors
   err = 1.0; ctr = 0 # This defines the change between the guess pin and it's neighbors and the count, these are used to determine convergence.
   while ((err > 1.0e-10)|(ctr<10000))&(ctr<200000):
      csp = (cs/(1+np.exp(-(1/500)*(5000-ctr))))+cs
      V[I_in+1,1] = -0.5
      V[Nx-I_out,-2] = 0.5  # 1
      V[0,1:Ny+1] = V[1,1:Ny+1]
      V[-1,1:Ny+1] = V[-2,1:Ny+1]
      V[1:Nx+1,0] = V[1:Nx+1,1]
      V[1:Nx+1,-1] = V[1:Nx+1,-2]
      
      dQ = csp*( (-2*V[1:Nx+1,1:Ny+1] + V[0:Nx,1:Ny+1] + V[2:Nx+2,1:Ny+1])/rxt \
         + (-2*V[1:Nx+1,1:Ny+1] + V[1:Nx+1,0:Ny] + V[1:Nx+1,2:Ny+2])/ryt  ) 
      
      V[1:Nx+1,1:Ny+1] += dQ
      err = np.abs(np.sum(dQ[1:Nx+1,1:Ny+1])/(Nx*Ny))
      ctr += 1
   Vlist = [T,rxt,ryt,dQ[I_in,0]/csp]
   for i in range(1,Nx+1):
      for j in range(1,Ny+1):
         Vlist.append(V[i,j])
   return Vlist


def convergex1(I_in, I_out, T, V):
   cs = convergence_speed_L(T) # Get convergence speed
   rxt = rx(T,Nx); ryt = newry(T,Ny) # Get strength of resistors
   err = 1.0; ctr = 0 # This defines the change between the guess pin and it's neighbors and the count, these are used to determine convergence.
   while ((err > 1.0e-12)|(ctr<10000))&(ctr<200000):
      csp = (cs/(1+np.exp(-(1/500)*(5000-ctr))))+cs
      V[I_in+1,1] = -0.5
      V[Nx-I_in,1] = 0.5  # 1
      V[0,1:Ny+1] = V[1,1:Ny+1]
      V[-1,1:Ny+1] = V[-2,1:Ny+1]
      V[1:Nx+1,0] = V[1:Nx+1,1]
      V[1:Nx+1,-1] = V[1:Nx+1,-2]
      
      dQ = csp*( (-2*V[1:Nx+1,1:Ny+1] + V[0:Nx,1:Ny+1] + V[2:Nx+2,1:Ny+1])/rxt \
         + (-2*V[1:Nx+1,1:Ny+1] + V[1:Nx+1,0:Ny] + V[1:Nx+1,2:Ny+2])/ryt  ) 
      
      V[1:Nx+1,1:Ny+1] += dQ
      err = np.abs(np.sum(dQ[1:Nx+1,1:Ny+1])/(Nx*Ny))
      ctr += 1
   Vlist = [T,rxt,ryt,dQ[I_in,0]/csp]
   for i in range(1,Nx+1):
      for j in range(1,Ny+1):
         Vlist.append(V[i,j])
   return Vlist


def convergex2(I_in, I_out, T, V):
   cs = convergence_speed_L(T) # Get convergence speed
   rxt = rx(T,Nx); ryt = newry(T,Ny) # Get strength of resistors
   err = 1.0; ctr = 0 # This defines the change between the guess pin and it's neighbors and the count, these are used to determine convergence.
   while ((err > 1.0e-12)|(ctr<10000))&(ctr<200000):
      csp = (cs/(1+np.exp(-(1/500)*(5000-ctr))))+cs
      V[I_out+1,-2] = -0.5
      V[Nx-I_out,-2] = 0.5  # 1
      V[0,1:Ny+1] = V[1,1:Ny+1]
      V[-1,1:Ny+1] = V[-2,1:Ny+1]
      V[1:Nx+1,0] = V[1:Nx+1,1]
      V[1:Nx+1,-1] = V[1:Nx+1,-2]
      
      dQ = csp*( (-2*V[1:Nx+1,1:Ny+1] + V[0:Nx,1:Ny+1] + V[2:Nx+2,1:Ny+1])/rxt \
         + (-2*V[1:Nx+1,1:Ny+1] + V[1:Nx+1,0:Ny] + V[1:Nx+1,2:Ny+2])/ryt  ) 
      
      V[1:Nx+1,1:Ny+1] += dQ
      err = np.abs(np.sum(dQ[1:Nx+1,1:Ny+1])/(Nx*Ny))
      ctr += 1
   Vlist = [T,rxt,ryt,dQ[I_out,-1]/csp]
   for i in range(1,Nx+1):
      for j in range(1,Ny+1):
         Vlist.append(V[i,j])
   return Vlist

def convergeL(T, V): #For new Geometry 27/03/23
   cs = 10 # Get convergence speed
   rxt = newry(T,Nx); ryt = rx(T,Ny) # Get strength of resistors
   err = 1.0; ctr = 0 # This defines the change between the guess pin and it's neighbors and the count, these are used to determine convergence.
   while ((err > 1.0e-20)|(ctr<10000))&(ctr<200000):
      csp = (cs/(1+np.exp(-(1/500)*(10000-ctr))))+cs
      V[1,2:-2] = -0.5 # Sets input and output pins to be constant
      V[-2,2:-2] = 0.5  # This used to be 1, and the other one used to be 1
      V[0,1:Ny+1] = V[1,1:Ny+1] # These four lines set the edge "ghost" pins to match their neighbors
      V[-1,1:Ny+1] = V[-2,1:Ny+1]
      V[1:Nx+1,0] = V[1:Nx+1,1]
      V[1:Nx+1,-1] = V[1:Nx+1,-2]
      
      # Calculating the change based on Ohm's law
      dQ = csp*( (-2*V[1:Nx+1,1:Ny+1] + V[0:Nx,1:Ny+1] + V[2:Nx+2,1:Ny+1])/rxt \
         + (-2*V[1:Nx+1,1:Ny+1] + V[1:Nx+1,0:Ny] + V[1:Nx+1,2:Ny+2])/ryt  ) 
      
      # Updating the matrix and convergence terms
      V[1:Nx+1,1:Ny+1] += dQ
      err = np.abs(np.sum(dQ[1:Nx+1,1:Ny+1])/(Nx*Ny))
      # print(err)
      ctr += 1
   Vlist = [T,rxt,ryt,dQ[0,1]/csp] # Setting the output lists
   for i in range(1,Nx+1):
      for j in range(1,Ny+1):
         Vlist.append(V[i,j])
   return Vlist

# Main simulate function
def simulate(I_in,I_out,type = 0):
   V = np.zeros([Nx+2,Ny+2]) # Empty Matrix
   data = [] # Output Data string
   N = 100 # Number of temperature steps.
   typestr = 'c' # typestr is the string of the type for output
   if type == 0:
      for Tctr in range(0,N):
         T = 300.0-Tctr*(300.-2.)/(N-1)
         data.append(convergec(I_in,I_out,T,V))
         V = np.zeros([Nx+2,Ny+2]) # Empty Matrix
   elif type == 1:
      typestr = 'd'
      for Tctr in range(0,N):
         T = 300.0-Tctr*(300.-2.)/(N-1)
         data.append(converged(I_in,I_out,T,V))
         V = np.zeros([Nx+2,Ny+2]) # Empty Matrix
   elif type == 2:
      typestr = 'x1'
      for Tctr in range(0,N):
         T = 300.0-Tctr*(300.-2.)/(N-1)
         data.append(convergex1(I_in,I_out,T,V))
         V = np.zeros([Nx+2,Ny+2]) # Empty Matrix
   elif type == 3:
      typestr = 'x2'
      for Tctr in range(0,N):
         T = 300.0-Tctr*(300.-2.)/(N-1)
         data.append(convergex2(I_in,I_out,T,V))
         V = np.zeros([Nx+2,Ny+2]) # Empty Matrix
   elif type == 4:
      typestr = 'L'
      for Tctr in range(0,N):
         T = 300.0-Tctr*(300.-2.)/(N-1)
         data.append(convergeL(T,V))
         V = np.zeros([Nx+2,Ny+2]) # Empty Matrix
   else:
      print('Error')

   headerlist = ['T','rx','ry','I'] # Header for the output file
   for i in range(1,Nx+1):
      for j in range(1,Ny+1):
         headerlist.append('V['+str(i)+','+str(j)+']')
   df = pd.DataFrame(data) # Making a dataframe from the data
   path = '../../Data/Sr327_Simulator/'+str(Nx)+'x'+str(Ny)+'DAT/' # Output path name
   if os.path.exists(path) == False:
      os.mkdir(path)
   df.to_csv(path+'Sr327_'+typestr+'_'+str(I_in)+'_to_'+str(I_out)+'.dat', index=False, header=headerlist)

   print('Simulation '+typestr+' '+str(I_in)+','+str(I_out)+' done')
   return True

if __name__ == "__main__":
   # # These lines define an arg-parser version that I was using previously, allowing you to use terminal to define variables.
   # parser = argparse.ArgumentParser(description='Plots 327 Simulation')
   # parser.add_argument("-c", "--save", help="Simulates all c-axis measurements", action="store_true")
   # parser.add_argument("-d", "--dimensions", help="Dimensions of simulation, x y = [x,y]", type=int, nargs=2, default=[12,3])
   # parser.add_argument("-p", "--pins", help="Guesses for input and output pins, i o = [i,o]", type=int, nargs=2, default=[3,4])
   # parser.add_argument("-e", "--experimental", help="Experimental version number", type=int, default=0)
   # args = parser.parse_args()

   # Nx, Ny = args.dimensions
   # guessinput, guessoutput = args.pins

   # # These lines define an sets of input and output pins.
   # input = []; output = []
   # for i in range(1, ceil(Nx/2)+1):
   #    input.append(i)
   #    output.append(i+1)
   # for i in range(2, ceil(Nx/2)+1):
   #    input.append(i)
   #    output.append(i-1)
   # for i in range(1, ceil(Nx/2)):
   #    input.append(i)
   #    output.append(i+2)
   # for i in range(3, ceil(Nx/2)+2):
   #    input.append(i)
   #    output.append(i-2)

   # # These process the simulations with multi-threading over the sets of input-output pins for a single type.
   # with Pool(processes=6) as pool:
   #    pool.starmap(simulate, product(input, repeat=2))
   # with Pool(processes=6) as pool:
   #    pool.starmap(simulate, zip(input,output,repeat(2)))
   # with Pool(processes=6) as pool2:
   #    pool2.starmap(simulate, zip(input,output,repeat(2)))

   # These process the simulations for a single guess of input-output pins, Pool is multi threading over types (c,diagonal, in plane).
   # guessinput = 5
   # guessoutput = 6
   # types = [0,1,2,3,4]
   # with Pool(processes=5) as pool:
   #   pool.starmap(simulate, zip(repeat(guessinput), repeat(guessoutput), types))
   # simulate(guessinput,guessoutput,2)

   # If you want to just simulate one particular set of guesses with a type use simulate() here
   simulate(5,6,type=4)