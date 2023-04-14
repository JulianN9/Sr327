
from email.mime import image
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from math import ceil
import argparse
import imageio.v2 as imageio
import os

fontsize = 24

def f(data,X,Y,T):
    Z = []
    z = []
    for y in Y:
        z = []
        for x in X:
            name = "V["+str(x)+","+str(y)+"]"
            z.append(data[name][T])
        Z.append(z)
    return Z

def RVTaxes(ax,ur = False):
    if ur == False:
        ax.set_xlabel(r'$T$  (K)',fontsize=fontsize)
        ax.set_ylabel(r'$\rho(T)$  (a.u.)',fontsize=fontsize)
        ax.legend(loc='upper right',prop={'size': 10})
        ax.set_xlim(0,300)
        ax.set_yticks([])
        ax.tick_params(axis='x', labelsize=16)
        ax.tick_params(axis='y', labelsize=16)
    else:
        ax.set_xlabel(r'$T$  (K)',fontsize=16)
        ax.set_ylabel(r'$\rho(T)$  (a.u.)',fontsize=16)
        ax.legend(loc='upper right')

def contouraxes(ax,heatmap,Nx,Ny,axtitle=True,gift = False):
    ax.set_xticks([1,ceil(Nx/2),Nx])
    ax.set_yticks([1,ceil(Ny/2),Ny])
    if gift == False:
        ax.tick_params(axis='x', labelsize=fontsize)
        ax.tick_params(axis='y', labelsize=fontsize)
        if axtitle == True:
            ax.set_xlabel('a-b plane',fontsize=fontsize)
            ax.set_ylabel('c-axis',fontsize=fontsize)
        cbr = plt.colorbar(heatmap,ax=ax)
        cbr.set_label('Voltage (V)', rotation=270,fontsize=fontsize)
        cbr.ax.get_yaxis().set_ticks([-0.5,0,0.5])
        cbr.ax.tick_params(labelsize=fontsize)
        cbr.ax.set_yticklabels(['-0.5','0','0.5'])
    else:
        ax.tick_params(axis='x')
        ax.tick_params(axis='y')
        if axtitle == True:
            ax.set_xlabel('a-b plane')
            ax.set_ylabel('c-axis')
        cbr = plt.colorbar(heatmap,ax=ax)
        cbr.set_label('Voltage (V)', rotation=270)
        cbr.ax.get_yaxis().set_ticks([-0.5,0,0.5])
        cbr.ax.set_yticklabels(['-0.5','0','0.5'])

def main(input,output,Nx,Ny,type = 'c',contours = 0,savefile = True, experimental = 0, gif = False):
    inputpin = "V["+str(input+1)+",1]"; outputpin = "V["+str(output+1)+","+str(Ny)+"]"
    mirrorinputpin = "V["+str(Nx-input)+",1]"; mirroroutputpin = "V["+str(Nx-output)+","+str(Ny)+"]"
    leftinput = "V[1,1]"; leftoutput = "V[1,"+str(Ny)+"]"
    middleinput = "V["+str(ceil(Nx/2))+",1]"; middleoutput = "V["+str(ceil(Nx/2))+","+str(Ny)+"]"
    righttinput = "V["+str(Nx)+",1]"; rightoutput = "V["+str(Nx)+","+str(Ny)+"]"
    files = []
    data = []
    path = '../../Data/Sr327_Simulator/'
    if experimental == 0:
        folder = path+str(Nx)+'x'+str(Ny)+'DAT/'       
    else:
        folder = path+'test'+str(experimental)+'_'+str(Nx)+'x'+str(Ny)+'DAT/'
    if (type == 'c') | (type == 'm'):
        files.append(folder+'Sr327_c_'+str(input)+'_to_'+str(output)+'.dat') 
    if (type == 'd') | (type == 'm'):
        files.append(folder+'Sr327_d_'+str(input)+'_to_'+str(output)+'.dat') 
    if (type == 'x1') | (type == 'm'):
        files.append(folder+'Sr327_x1_'+str(input)+'_to_'+str(output)+'.dat') 
    if (type == 'x2'):
        files.append(folder+'Sr327_x2_'+str(input)+'_to_'+str(output)+'.dat') 
    if (type == 'L'):
        files.append(folder+'Sr327_L_'+str(input)+'_to_'+str(output)+'.dat') 
    for fname in files:
        data.append(pd.read_csv(fname,header=[0])) 

    if (contours == 0)&(gif==False):
        fig = plt.figure(figsize=(8,6))
        ax = fig.add_axes([0.2,0.2,0.7,0.6])

        if (type =='m'):
            ax.scatter(data[0]["T"],data[0]["rx"],color='red',s=1.0,label='rx(T)')
            # ax.scatter(data[0]["T"],0.4*data[0]["ry"],color='blue',s=1.0,label='ry(T) [0.4]')
            ax.scatter(data[0]["T"],data[0]["ry"],color='blue',s=1.0,label='ry(T)')
        if (type == 'c')| (type == 'd'):
            ax.scatter(data[0]["T"],2.65*(data[0][leftoutput]-data[0][leftinput])/data[0]["I"],color='purple',s=8.0,marker='.',label='Left Edge') 
            ax.scatter(data[0]["T"],2.65*(data[0][middleoutput]-data[0][middleinput])/data[0]["I"],color='green',s=8.0,marker='+',label='Center') 
            ax.scatter(data[0]["T"],2.65*(data[0][rightoutput]-data[0][righttinput])/data[0]["I"],color='brown',s=8.0,marker='x',label='Right Edge')
        if type == 'c':
            ax.scatter(data[0]["T"],2.65*(data[0][outputpin]-data[0][inputpin])/data[0]["I"],color='black',s=8.0,marker='v',label='Input')
            ax.scatter(data[0]["T"],2.65*(data[0][mirroroutputpin]-data[0][mirrorinputpin])/data[0]["I"],color='yellow',s=8.0,marker='^',label='Output')
            # ax.scatter(data[0]["T"],2.65*(data[0][mirroroutputpin]-data[0][mirrorinputpin]),color='yellow',s=8.0,marker='^',label='Voltage')
            # ax.scatter(data[0]["T"],data[0]["I"],color='purple',s=8.0,marker='.',label='Current') 
        if type == 'd':
            ax.scatter(data[0]["T"],2.65*(data[0][mirroroutputpin ]-data[0][inputpin])/data[0]["I"],color='black',s=8.0,marker='v',label='Input')
            ax.scatter(data[0]["T"],2.65*(data[0][outputpin]-data[0][mirrorinputpin])/data[0]["I"],color='yellow',s=8.0,marker='^',label='Output')
        if (type == 'x1')|(type == 'x2'):
            ax.scatter(data[0]["T"],2.65*(data[0][rightoutput]-data[0][leftoutput])/data[0]["I"],color='purple',s=8.0,marker='.',label='Top Edge') 
            ax.scatter(data[0]["T"],2.65*(data[0][righttinput]-data[0][leftinput])/data[0]["I"],color='brown',s=8.0,marker='x',label='Bottom Edge')
            ax.scatter(data[0]["T"],2.65*(data[0][mirrorinputpin]-data[0][inputpin])/data[0]["I"],color='black',s=8.0,marker='v',label='Input')
            ax.scatter(data[0]["T"],2.65*(data[0][mirroroutputpin]-data[0][outputpin])/data[0]["I"],color='yellow',s=8.0,marker='^',label='Output')
        if type == 'L':
            # ax.scatter(data[0]["T"],data[0]["I"],s=8.0,marker='.',label='Current') 
            # for x in range(1,8):
            #     ax.scatter(data[0]["T"],-2.65*(data[0]["V["+str(1+x)+",1]"]-data[0]["V["+str(Nx-x)+",1]"])/data[0]["I"],s=8.0,marker='.',label='Top Edge '+str(x)) 
            #     ax.scatter(data[0]["T"],2.65*(data[0]["V["+str(1+x)+",1]"]-data[0]["V["+str(Nx-x)+",1]"]),s=8.0,marker='.',label='Voltage '+str(x)) 
            ax.scatter(data[0]["T"],-2.65*(data[0]["V[8,1]"]-data[0]["V["+str(Nx-7)+",1]"])/data[0]["I"],color='purple',s=8.0,marker='.',label='Top Edge') 
        if type == 'm':
            ax.plot(data[0]["T"],2.65*(data[0][mirroroutputpin]-data[0][mirrorinputpin])/data[0]["I"],color='orange',marker='^',label='c-axis')
            ax.plot(data[1]["T"],2.65*(data[1][outputpin]-data[1][mirrorinputpin])/data[1]["I"],color='green',marker='^',label='diagonal')
            ax.plot(data[2]["T"],2.65*(data[2][mirroroutputpin]-data[2][outputpin])/data[2]["I"],color='red',marker='^',label='in-plane')
            #ax.plot(data[2]["T"],2.65*(data[2][mirrorinputpin]-data[2][inputpin])/data[2]["I"],color='red',marker='^')
        RVTaxes(ax)
        plt.show()
        if savefile == True:
            if experimental > 0:
                fig.savefig('../../Plots/Sr327/Simulations/test'+str(experimental)+'_'+str(Nx)+'x'+str(Ny)+'_'+type+'.svg')
            else:
                fig.savefig('../../Plots/Sr327/Simulations/'+str(Nx)+'x'+str(Ny)+'_'+type+'.svg')

    elif gif == False:
        x = []; y = []
        for i in range(1,Nx+1):
            x.append(i)
        for j in range(1,Ny+1):
            y.append(j)
        X, Y = np.meshgrid(x,y)

        Z = f(data[0],x,y,contours-1)
        fig = plt.figure(figsize=(6,4))
        ax = fig.add_axes([0.2,0.2,0.7,0.6])
        heatmap = ax.contourf(X,Y, Z,cmap='RdGy')

        contouraxes(ax,heatmap,Nx,Ny)

        plt.show()

        if (savefile == True)&(gif == False):
            if experimental > 0:
                fig.savefig('../../Plots/Sr327/Simulations/test'+str(experimental)+'_'+str(Nx)+'x'+str(Ny)+'_'+type+'_contour_'+str(contours)+'.svg')
            else:
                fig.savefig('../../Plots/Sr327/Simulations/'+str(Nx)+'x'+str(Ny)+'_'+type+'_contour_'+str(contours)+'.svg')

    else:
        x = []; y = []
        for i in range(1,Nx+1):
            x.append(i)
        for j in range(1,Ny+1):
            y.append(j)
        X, Y = np.meshgrid(x,y)
        gifdir = 'tempgif/'
        os.mkdir(gifdir)
        images = []
        images2 = []
        if not(type == 'm'):
            for i in range(0,100):
                Z = f(data[0],x,y,i)
                fig = plt.figure(figsize=(6,4))
                ax = fig.add_axes([0.2,0.2,0.7,0.6])
                heatmap = ax.contourf(X,Y, Z,cmap='RdGy')
                contouraxes(ax,heatmap,Nx,Ny,gift=True)

                fig.savefig(gifdir+str(i)+".png")
                images.append(imageio.imread(gifdir+str(i)+".png"))     
                plt.close(fig=fig)       
        else:
            for i in range(0,100):
                fig = plt.figure(figsize=(8,5))
                grid = plt.GridSpec(6,4, hspace=0.2, wspace=0.2)
                main_ax = fig.add_subplot(grid[:-1,:-1])
            
                c_ax = fig.add_subplot(grid[0:1,-1:])
                d_ax = fig.add_subplot(grid[2:-3,-1:])
                x_ax = fig.add_subplot(grid[-2:-1,-1:])
                contour_axes = [c_ax,d_ax,x_ax]
                for number, ax in enumerate(contour_axes):
                    Z = f(data[number],x,y,i)
                    if number+1 == 1:
                        ax.set_title('c-axis')
                    if number+1 == 2:
                        ax.set_title('diagonal')
                    if number+1 == 3:
                        ax.set_title('in-plane')
                    if number+1 == 4:
                        ax.set_title('long c-axis')
                    heatmap = ax.contourf(X,Y, Z,cmap='RdGy')
                    contouraxes(ax,heatmap,Nx,Ny,False,True)

                main_ax.plot(data[0]["T"][:i+1],2.65*(data[0][mirroroutputpin][:i+1]-data[0][mirrorinputpin][:i+1])/data[0]["I"][:i+1],color='orange',marker='^',label='c-axis')
                main_ax.plot(data[1]["T"][:i+1],2.65*(data[1][outputpin][:i+1]-data[1][mirrorinputpin][:i+1])/data[1]["I"][:i+1],color='green',marker='^',label='diagonal')
                main_ax.plot(data[2]["T"][:i+1],2.65*(data[2][mirroroutputpin][:i+1]-data[2][outputpin][:i+1])/data[2]["I"][:i+1],color='red',marker='^',label='in-plane')
                RVTaxes(main_ax,True)
                fig.savefig(gifdir+str(i)+".png")
                images.append(imageio.imread(gifdir+str(i)+".png"))
                plt.close(fig=fig)       
            ticker = 0
            while ticker < 50:
                images.append(imageio.imread(gifdir+str(99)+".png"))
                ticker += 1
        imageio.mimsave('../../Plots/Sr327/Simulations/gif_'+str(Nx)+'x'+str(Ny)+'_'+type+'.gif',images,fps=10)
        for i in range(0,100):
            os.remove(gifdir+str(i)+".png")
        os.rmdir(gifdir)
    return True

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Plots 327 Simulation')
    parser.add_argument("-t", "--type", help="Which type, c = c-axis, d = diagonal, x1,2 = in-plane or m = mixed")
    parser.add_argument("-c", "--contours", help="Plot contours instead of resistivity, value is temperature", type=int, default=0)
    parser.add_argument("-s", "--save", help="Enables saving the plot", action="store_true")
    parser.add_argument("-d", "--dimensions", help="Dimensions of simulation, x y = [x,y]", type=int, nargs=2, default=[24,5])
    parser.add_argument("-p", "--pins", help="Input and output pins, i o = [i,o]", type=int, nargs=2, default=[2,3])
    parser.add_argument("-e", "--experimental", help="Experimental version number", type=int, default=0)
    parser.add_argument("-g", "--gif", help="Enables giffing the plot", action="store_true")
    args = parser.parse_args()
    if (args.type == 'c')|(args.type == 'd')|(args.type == 'm')|(args.type == 'x1')|(args.type == 'x2')|(args.type == 'L'):
        main(args.pins[0],args.pins[1],args.dimensions[0],args.dimensions[1],type=args.type,contours=args.contours,savefile=args.save,experimental=args.experimental,gif=args.gif)
    else:
        print("No type specified, defaulting to c-axis")
        main(args.pins[0],args.pins[1],args.dimensions[0],args.dimensions[1],contours=args.contours,savefile=args.save,experimental=args.experimental,gif=args.gif)
    