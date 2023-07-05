
from email.mime import image
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from math import ceil
from pathlib import Path
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
    ax.set_xlabel(r'Temperature  (K)',fontsize=fontsize)
    ax.set_ylabel(r'$\rho(T)$  (a.u.)',fontsize=fontsize)
    ax.set_yticks([])
    ax.legend(loc='upper right',prop={'size': 10})
    ax.tick_params(axis='x', labelsize=16)
    if ur == False:
        ax.set_xlim(0,300)
        ax.legend(fontsize='16')

def contouraxes(ax,heatmap,Nx,Ny,axtitle=1,gift = False,mixed = 'none'):
    ax.set_xticks([1,ceil(Nx/2),Nx])
    ax.set_yticks([1,ceil(Ny/2),Ny])
    if (gift == False)&(mixed == 'none'):
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
    # elif(gift==True):
    #     ax.set_xticks([])
    #     ax.set_yticks([])
    #     ax.tick_params(axis='x')
    #     ax.tick_params(axis='y')
    #     if axtitle > 0:
    #         ax.set_xlabel('a-b plane')
    #         ax.set_ylabel('c-axis')
        # cbr = plt.colorbar(heatmap,ax=ax)
        # cbr.set_label('Voltage (V)', rotation=270)
        # cbr.ax.get_yaxis().set_ticks([-0.5,0,0.5])
        # cbr.ax.set_yticklabels(['-0.5','0','0.5'])
    else:
        if axtitle == 2:
            ax.set_ylabel('a-b plane')
            ax.set_xlabel('c-axis')
        elif axtitle > 0:
            ax.set_xlabel('a-b plane')
            ax.set_ylabel('c-axis')
        if(mixed == 'top'):
            ax.set_xticks([])
            ax.set_yticks([])
        if(mixed == 'left'):
            ax.set_xticks([])
            if (gift==False):
                ax.set_ylabel(r'$\bf{300 K}$' '\n' 'c-axis')
        if(mixed == 'corner'):
            if (gift==False):
                ax.set_ylabel(r'$\bf{2 K}$' '\n' 'c-axis')
        if(mixed == 'bottom'):
            ax.set_yticks([])

def loadsimdata(type,Nx,Ny,input,output,experimental = 0,pressure = 0,implicit = False):
    files = []
    data = []
    path = '../../Data/Sr327_Simulator/'
    if implicit == True:
        path = '../../Data/Sr327_ImplicitSimulator/' 
    if experimental == 0:
        if pressure == 0:
            folder = path+str(Nx)+'x'+str(Ny)+'DAT/'       
        else:
            folder = path+str(pressure)+'P'+str(Nx)+'x'+str(Ny)+'DAT/'
    else:
        folder = path+'test'+str(experimental)+'_'+str(Nx)+'x'+str(Ny)+'DAT/'
    L_check = Path(folder+'Sr327_L_'+str(input)+'_to_'+str(output)+'.dat')
    if (type == 'c') | (type == 'm'):
        files.append(folder+'Sr327_c_'+str(input)+'_to_'+str(output)+'.dat') 
    if (type == 'd') | (type == 'm'):
        files.append(folder+'Sr327_d_'+str(input)+'_to_'+str(output)+'.dat') 
    if (type == 'x1') | (type == 'm'):
        files.append(folder+'Sr327_x1_'+str(input)+'_to_'+str(output)+'.dat') 
    if (type == 'x2'):
        files.append(folder+'Sr327_x2_'+str(input)+'_to_'+str(output)+'.dat') 
    if (type == 'L')| ((type == 'm')&(L_check.is_file() == True)):
        files.append(folder+'Sr327_L_'+str(input)+'_to_'+str(output)+'.dat') 
    for fname in files:
        data.append(pd.read_csv(fname,header=[0])) 
    return data, L_check

def main(input,output,Nx,Ny,type = 'c',contours = 0,savefile = True, experimental = 0, gif = False, pressure = 0, implicit = False):
    inputpin = "V["+str(input+1)+",1]"; outputpin = "V["+str(output+1)+","+str(Ny)+"]"
    mirrorinputpin = "V["+str(Nx-input)+",1]"; mirroroutputpin = "V["+str(Nx-output)+","+str(Ny)+"]"
    leftinput = "V[1,1]"; leftoutput = "V[1,"+str(Ny)+"]"
    middleinput = "V["+str(ceil(Nx/2))+",1]"; middleoutput = "V["+str(ceil(Nx/2))+","+str(Ny)+"]"
    righttinput = "V["+str(Nx)+",1]"; rightoutput = "V["+str(Nx)+","+str(Ny)+"]"
    data, L_check = loadsimdata(type,Nx,Ny,input,output,experimental,pressure,implicit)

    if (contours == 0)&(gif==False):
        fig = plt.figure(figsize=(10,8))
        ax = fig.add_axes([0.1,0.1,0.85,0.85])

        if (type =='m')|(type == 'c'):
            ax.scatter(data[0]["T"],0.36*data[0]["rx"],color='red',s=4.0,label= r'$\rho_{ab}$(T)')
            ax.scatter(data[0]["T"],0.36*data[0]["ry"],color='blue',s=4.0,label= r'$\rho_c$(T)')
            # ax.scatter(data[0]["T"],data[0]["ry"],color='blue',s=1.0,label='ry(T)')
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
            ax.scatter(data[0]["T"],2.65*(data[0][mirroroutputpin]-data[0][outputpin])/data[0]["I"],color='red',s=8.0,marker='^',label='Output')
            # ax.scatter(data[0]["T"],data[0]["rx"],color='red',s=4.0,label= r'$\rho_{ab}$(T)')
        if type == 'L':
            # ax.scatter(data[0]["T"],data[0]["I"],s=8.0,marker='.',label='Current') 
            for x in range(1,8):
                ax.scatter(data[0]["T"],-2.65*(data[0]["V["+str(1+x)+",1]"]-data[0]["V["+str(Nx-x)+",1]"])/data[0]["I"],s=8.0,marker='.',label='long c-axis '+str(x)) 
                # ax.scatter(data[0]["T"],2.65*(data[0]["V["+str(1+x)+",1]"]-data[0]["V["+str(Nx-x)+",1]"]),s=8.0,marker='.',label='Voltage '+str(x)) 
            # ax.scatter(data[0]["T"],-2.65*(data[0]["V[8,1]"]-data[0]["V["+str(Nx-7)+",1]"])/data[0]["I"],color='purple',s=8.0,marker='.',label='long c-axis') 
        if type == 'm':
            scale = 0.36*data[0]["ry"][99]
            ax.plot(data[0]["T"],(scale/((data[0][mirroroutputpin][99]-data[0][mirrorinputpin][99])/data[0]["I"][99]))*(data[0][mirroroutputpin]-data[0][mirrorinputpin])/data[0]["I"],color='C1',marker='^',label='c-axis')
            ax.plot(data[1]["T"],(scale/((data[1][outputpin][99]-data[1][mirrorinputpin][99])/data[1]["I"][99]))*(data[1][outputpin]-data[1][mirrorinputpin])/data[1]["I"],color='green',marker='^',label='diagonal')
            ax.plot(data[2]["T"],(Ny/(Nx-output-input))*2.65*(data[2][mirroroutputpin]-data[2][outputpin])/data[2]["I"],color='red',marker='^',label='in-plane')
            #ax.plot(data[2]["T"],2.65*(data[2][mirrorinputpin]-data[2][inputpin])/data[2]["I"],color='red',marker='^')
            if L_check.is_file() == True:
                ax.plot(data[3]["T"],(scale/((data[3]["V[8,1]"][99]-data[3]["V["+str(Nx-7)+",1]"][99])/data[3]["I"][99]))*(data[3]["V[8,1]"]-data[3]["V["+str(Nx-7)+",1]"])/data[3]["I"],color='purple',marker='^',label='long c-axis') 
        RVTaxes(ax)
        plt.show()
        if (savefile == True):
            if implicit == False:
                if experimental > 0:
                    fig.savefig('../../Plots/Sr327/Simulations/test'+str(experimental)+'_'+str(Nx)+'x'+str(Ny)+'_'+type+'.svg')
                else:
                    fig.savefig('../../Plots/Sr327/Simulations/'+str(Nx)+'x'+str(Ny)+'_'+type+'.svg')
            else:
                if pressure != 0:
                    fig.savefig('../../Plots/Sr327/ImplicitSimulations/'+str(pressure)+'P'+str(Nx)+'x'+str(Ny)+'_'+type+'.svg')
                else:
                    fig.savefig('../../Plots/Sr327/ImplicitSimulations/'+str(Nx)+'x'+str(Ny)+'_'+type+'.svg')

    elif gif == False:
        x = []; y = []
        for i in range(1,Nx+1):
            x.append(i)
        for j in range(1,Ny+1):
            y.append(j)
        X, Y = np.meshgrid(x,y)
        if (type == 'm')&(L_check.is_file() == True):
            fig = plt.figure(figsize=(16,5))
            grid = plt.GridSpec(4,16, hspace=0.4, wspace=0.8)
            c_low_ax = fig.add_subplot(grid[:-2,0:4])
            d_low_ax = fig.add_subplot(grid[:-2,4:8])
            x_low_ax = fig.add_subplot(grid[:-2,8:12])
            L_low_ax = fig.add_subplot(grid[:-2,12:])
            c_high_ax = fig.add_subplot(grid[-2:,0:4])
            d_high_ax = fig.add_subplot(grid[-2:,4:8])
            x_high_ax = fig.add_subplot(grid[-2:,8:12])
            L_high_ax = fig.add_subplot(grid[-2:,12:])
            contour_axes = [[c_low_ax,c_high_ax],[d_low_ax,d_high_ax],[x_low_ax,x_high_ax],[L_low_ax,L_high_ax]]
            for number, ax in enumerate(contour_axes):
                heatmap1 = ax[0].contourf(X,Y, f(data[number],x,y,0),cmap='RdGy')
                heatmap2 = ax[1].contourf(X,Y, f(data[number],x,y,99),cmap='RdGy')
                if number+1 == 1:
                    ax[0].set_title('c-axis',fontweight='bold')
                    contouraxes(ax[0],heatmap1,Nx,Ny,1,mixed='left')
                    contouraxes(ax[1],heatmap2,Nx,Ny,1,mixed='corner')
                    cbr = fig.colorbar(heatmap1,ax=contour_axes)
                    cbr.set_label('Voltage (V)', rotation=270,fontweight='bold')
                    cbr.ax.get_yaxis().set_ticks([-0.5,0,0.5])
                    cbr.ax.set_yticklabels(['-0.5','0','0.5'])
                if number+1 == 2:
                    ax[0].set_title('diagonal',fontweight='bold')
                    contouraxes(ax[0],heatmap1,Nx,Ny,1,mixed='top')
                    contouraxes(ax[1],heatmap2,Nx,Ny,1,mixed='bottom')
                if number+1 == 3:
                    ax[0].set_title('in-plane',fontweight='bold')
                    contouraxes(ax[0],heatmap1,Nx,Ny,1,mixed='top')
                    contouraxes(ax[1],heatmap2,Nx,Ny,1,mixed='bottom')
                if number+1 == 4:
                    ax[0].set_title('long c-axis',fontweight='bold')
                    contouraxes(ax[0],heatmap1,Nx,Ny,2,mixed='top')
                    contouraxes(ax[1],heatmap2,Nx,Ny,2,mixed='bottom')
            plt.show()
        else:
            Z = f(data[0],x,y,contours-1)
            fig = plt.figure(figsize=(6,4))
            ax = fig.add_axes([0.2,0.2,0.7,0.6])
            heatmap = ax.contourf(X,Y, Z,cmap='RdGy')

            contouraxes(ax,heatmap,Nx,Ny)

            plt.show()

        if (savefile == True):
            if implicit == False:
                if experimental > 0:
                    fig.savefig('../../Plots/Sr327/Simulations/test'+str(experimental)+'_'+str(Nx)+'x'+str(Ny)+'_'+type+'_contour_'+str(contours)+'.svg')
                else:
                    fig.savefig('../../Plots/Sr327/Simulations/'+str(Nx)+'x'+str(Ny)+'_'+type+'_contour_'+str(contours)+'.svg')
            else:
                if pressure != 0:
                    fig.savefig('../../Plots/Sr327/ImplicitSimulations/'+str(pressure)+'P'+str(Nx)+'x'+str(Ny)+'_'+type+'_contour_'+str(contours)+'.svg')
                else:
                    fig.savefig('../../Plots/Sr327/ImplicitSimulations/'+str(Nx)+'x'+str(Ny)+'_'+type+'_contour_'+str(contours)+'.svg')


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
            if L_check.is_file() == True:
                for i in range(0,100):
                    fig = plt.figure(figsize=(12,6))
                    grid = plt.GridSpec(8,4, hspace=0.2, wspace=0.2)
                    main_ax = fig.add_subplot(grid[:-1,:-1])
                
                    c_ax = fig.add_subplot(grid[0:1,-1:])
                    d_ax = fig.add_subplot(grid[2:3,-1:])
                    x_ax = fig.add_subplot(grid[4:-3,-1:])
                    L_ax = fig.add_subplot(grid[-2:-1,-1:])
                    contour_axes = [c_ax,d_ax,x_ax,L_ax]
                    for number, ax in enumerate(contour_axes):
                        Z = f(data[number],x,y,i)
                        heatmap = ax.contourf(X,Y, Z,cmap='RdGy')
                        if number+1 == 1:
                            ax.set_title('c-axis')
                            cbr = fig.colorbar(heatmap,ax=contour_axes)
                            cbr.set_label('Voltage (V)', rotation=270,fontweight='bold')
                            cbr.ax.get_yaxis().set_ticks([-0.5,0,0.5])
                            cbr.ax.set_yticklabels(['-0.5','0','0.5'])
                            contouraxes(ax,heatmap,Nx,Ny,1,mixed='top')
                        if number+1 == 2:
                            ax.set_title('diagonal')
                            contouraxes(ax,heatmap,Nx,Ny,1,mixed='top')
                        if number+1 == 3:
                            ax.set_title('in-plane')
                            contouraxes(ax,heatmap,Nx,Ny,1,mixed='top')
                        if number+1 == 4:
                            ax.set_title('long c-axis')
                            contouraxes(ax,heatmap,Nx,Ny,2,mixed='top')
                        # contouraxes(ax,heatmap,Nx,Ny,0,True)
                        
                    scale = 0.36*data[0]["ry"][99]
                    main_ax.plot(data[0]["T"][:i+1],(scale/((data[0][mirroroutputpin][99]-data[0][mirrorinputpin][99])/data[0]["I"][99]))*(data[0][mirroroutputpin][:i+1]-data[0][mirrorinputpin][:i+1])/data[0]["I"][:i+1],color='orange',marker='^',label='c-axis')
                    main_ax.plot(data[1]["T"][:i+1],(scale/((data[1][outputpin][99]-data[1][mirrorinputpin][99])/data[1]["I"][99]))*(data[1][outputpin][:i+1]-data[1][mirrorinputpin][:i+1])/data[1]["I"][:i+1],color='green',marker='^',label='diagonal')
                    main_ax.plot(data[2]["T"][:i+1],2.65*(data[2][mirroroutputpin][:i+1]-data[2][outputpin][:i+1])/data[2]["I"][:i+1],color='red',marker='^',label='in-plane')
                    main_ax.plot(data[3]["T"][:i+1],(scale/((data[3]["V[8,1]"][99]-data[3]["V["+str(Nx-7)+",1]"][99])/data[3]["I"][99]))*(data[3]["V[8,1]"][:i+1]-data[3]["V["+str(Nx-7)+",1]"][:i+1])/data[3]["I"][:i+1],color='purple',marker='^',label='long c-axis') 
                    RVTaxes(main_ax,True)
                    fig.savefig(gifdir+str(i)+".png")
                    images.append(imageio.imread(gifdir+str(i)+".png"))
                    plt.close(fig=fig)       
                ticker = 0
                while ticker < 50:
                    images.append(imageio.imread(gifdir+str(99)+".png"))
                    ticker += 1
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
                        heatmap = ax.contourf(X,Y, Z,cmap='RdGy')
                        contouraxes(ax,heatmap,Nx,Ny,0,True)

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
        if implicit == False:
            imageio.mimsave('../../Plots/Sr327/Simulations/gif_'+str(Nx)+'x'+str(Ny)+'_'+type+'.gif',images,fps=10)
        else:
            if pressure == 0:
                imageio.mimsave('../../Plots/Sr327/ImplicitSimulations/gif_'+str(Nx)+'x'+str(Ny)+'_'+type+'.gif',images,fps=10)
            else:
                imageio.mimsave('../../Plots/Sr327/ImplicitSimulations/gif_'+str(pressure)+'P'+str(Nx)+'x'+str(Ny)+'_'+type+'.gif',images,fps=10)
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
    parser.add_argument("-P", "--pressure", help="Selects which pressure to use for implicit sim", type=int, default=0)
    args = parser.parse_args()
    if (args.type == 'c')|(args.type == 'd')|(args.type == 'm')|(args.type == 'x1')|(args.type == 'x2')|(args.type == 'L'):
        main(args.pins[0],args.pins[1],args.dimensions[0],args.dimensions[1],type=args.type,contours=args.contours,savefile=args.save,experimental=args.experimental,gif=args.gif,pressure=args.pressure,implicit=True)
    else:
        print("No type specified, defaulting to c-axis")
        main(args.pins[0],args.pins[1],args.dimensions[0],args.dimensions[1],contours=args.contours,savefile=args.save,experimental=args.experimental,gif=args.gif)
    