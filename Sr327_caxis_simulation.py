import argparse
import os
import numpy as np
import pandas as pd
from math import exp, ceil
from multiprocessing import Pool
from itertools import product, repeat

#You should be using simulation.py

def R_c(T):  # in mu-Ohm-cm (hence factor of 1000.0)
  f = 1.0/(exp((T-10.0)/10.) + 1.)
  g = 1.0/(exp((T-40.0)/20.) + 1.)
  return 1000.0*(f*(1.0+0.02*T*T) + 0.08*T*g*(1.-f) + (8.0+0.0005*T)*(1.-g)*(1.-f))

def R_ab(T): 
  f = 1.0/(exp((T-20.)/10.) + 1.0)
  return (1.7+0.03*T*T)*f + 0.68*T*(1.0-f)

rx = 2.30 # x direction resistor
ry = 2000.0 # y direction resistor

Nx = 24 # x nodes
Ny = 6 # y nodes
# I_in = 6 
# I_out = 7  # Nx-1 
# type = 'd'

def simulate(I_in, I_out,typeint = 0):
   type = 'c'
   V = np.zeros([Nx+2,Ny+2])
   if typeint == 0:
      type = 'c'
   if typeint == 1:
      type = 'd'
   if typeint == 2:
      type = 'x1'
   if typeint == 3:
      type = 'x2'

   test=0
   contours=0
   ## for j in range(0,20):    # range(19,20)  # for contour plot
   ##    rx = 2.0 + 1.0*j
   data = []
   N = 100
   for Tctr in range(0,N): # 50):
      T = 300.0-Tctr*(300.-2.)/(N-1)
      convergence_speed = 10
      if T < 50:
         convergence_speed = 1
      if T < 10:
         convergence_speed = 0.5
      # rx = 20*R_ab(T)/(Nx-2)
      # ry = 20*R_c(T)/(Ny-2)
      # rx = 440*R_ab(T)/(Nx-1)
      # ry = 20*R_c(T)/(Ny-1)
      rx = 200*R_ab(T)/(Nx-1)
      ry = 20*R_c(T)/(Ny-1)
      # rx = 2*R_ab(T)/(Nx-1)
      # ry = R_c(T)/(Ny-1)
      delta_Q = 1.0; ctr = 0
      # for i in range(0,2000000):
      while ((delta_Q > 1.0e-6)|(ctr<10000))&(ctr<1000000):
         if type == 'c':
            V[I_in+1,1] = -0.5 
            V[I_out+1,-2] = 0.5  # 1
         if type == 'd':
            V[I_in+1,1] = -0.5
            V[Nx-I_out,-2] = 0.5  # 1
         if type == 'x1':
            V[I_in+1,1] = -0.5
            V[Nx-I_in,1] = 0.5  # 1
         if type == 'x2':
            V[I_out+1,-2] = -0.5
            V[Nx-I_out,-2] = 0.5  # 1
         V[0,1:Ny+1] = V[1,1:Ny+1]
         V[-1,1:Ny+1] = V[-2,1:Ny+1]
         V[1:Nx+1,0] = V[1:Nx+1,1]
         V[1:Nx+1,-1] = V[1:Nx+1,-2]
      
         if ctr%9900 == 0:
            if test==1:
               for j in range(0,Nx+2):
                  print ("{:2.3f} {:2.3f} {:2.3f} {:2.3f} {:2.3f} {:2.3f} {:2.3f} {:2.3f} {:2.3f} {:2.3f} {:2.3f} {:2.3f}".format(V[j,0],V[j,1],V[j,2],V[j,3],V[j,4],V[j,5],V[j,6],V[j,7],V[j,8],V[j,9],V[j,10],V[j,11]))
               print (' ')
         
         dQ = convergence_speed*( (-2*V[1:Nx+1,1:Ny+1] + V[0:Nx,1:Ny+1] + V[2:Nx+2,1:Ny+1])/rx \
            + (-2*V[1:Nx+1,1:Ny+1] + V[1:Nx+1,0:Ny] + V[1:Nx+1,2:Ny+2])/ry  ) 
         
         V[1:Nx+1,1:Ny+1] += dQ
         if type == 'c':
            delta_Q = (dQ[I_in,0] + dQ[I_out,-1])
         if type == 'd':
            delta_Q = (dQ[I_in,0] + dQ[(Nx+1-(I_out)),-1])
         if type == 'x1':
            delta_Q = (dQ[I_in,0] + dQ[(Nx+1-(I_in)),0])
         if type == 'x2':
            delta_Q = (dQ[I_out,-1] + dQ[(Nx+1-(I_out)),-1])
         # if ctr%20000==0:
         #   print ('#',ctr,delta_Q,dQ[I_in,0],dQ[I_out,-1])
         ctr += 1
      if contours==0:
         Vlist = [T,rx,ry,dQ[I_in,0]/convergence_speed]
         for i in range(1,Nx+1):
            for j in range(1,Ny+1):
               Vlist.append(V[i,j])
         data.append(Vlist)

   if contours == 0:
      headerlist = ['T','rx','ry','I']
      for i in range(1,Nx+1):
         for j in range(1,Ny+1):
            headerlist.append('V['+str(i)+','+str(j)+']')
      df = pd.DataFrame(data)
      path = 'test11_'+str(Nx)+'x'+str(Ny)+'DAT/'
      if os.path.exists(path) == False:
         os.mkdir(path)
      df.to_csv(path+'Sr327_'+type+'_'+str(I_in)+'_to_'+str(I_out)+'.dat', index=False, header=headerlist)

   if contours==1:
      for i in range(1,Nx+1):
         for j in range(1,Ny+1):
            print (i,j,V[i,j])
         print (' ')

   print('Simulation '+type+' '+str(I_in)+','+str(I_out)+' done')
   return 0

if __name__ == "__main__":
   # parser = argparse.ArgumentParser(description='Plots 327 Simulation')
   # parser.add_argument("-c", "--save", help="Simulates all c-axis measurements", action="store_true")
   # parser.add_argument("-d", "--dimensions", help="Dimensions of simulation, x y = [x,y]", type=int, nargs=2, default=[12,3])
   # parser.add_argument("-p", "--pins", help="Guesses for input and output pins, i o = [i,o]", type=int, nargs=2, default=[3,4])
   # parser.add_argument("-e", "--experimental", help="Experimental version number", type=int, default=0)
   # args = parser.parse_args()

   # Nx, Ny = args.dimensions
   # guessinput, guessoutput = args.pins

   input = []; output = []
   for i in range(1, ceil(Nx/2)+1):
      input.append(i)
      output.append(i+1)
   # for i in range(2, ceil(Nx/2)+1):
   #    input.append(i)
   #    output.append(i-1)
   # for i in range(1, ceil(Nx/2)):
   #    input.append(i)
   #    output.append(i+2)
   # for i in range(3, ceil(Nx/2)+2):
   #    input.append(i)
   #    output.append(i-2)
   # with Pool(processes=6) as pool:
   #    pool.starmap(simulate, product(input, repeat=2))
   with Pool(processes=6) as pool:
      pool.starmap(simulate, zip(input,output,repeat(2)))
   # with Pool(processes=6) as pool2:
   #    pool2.starmap(simulate, zip(input,output,repeat(2)))
   # guessinput = 7
   # guessoutput = 8
   # types = [0,1,2,3]
   # with Pool(processes=4) as pool:
   #   pool.starmap(simulate, zip(repeat(guessinput), repeat(guessoutput), types))
   #simulate(guessinput,guessoutput,2)