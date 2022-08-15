"""Analysis file of 3d rotating rayleigh benard convection

Usage:
    plot-txt-file.py --dir=<directory> [--t=<transient> --fig=<Figure>]
    plot-txt-file.py -h | --help

Options:
    -h --help   		          	 Display this help message
    --dir=<directory>        		 Directory
    --t=<transient>   			     Transient to be ignored [default: 2000]
    --fig=<Figure>    			     Produce Figures [default: 1]
"""

from docopt import docopt
import numpy as np
import matplotlib.pyplot as plt
import os
import h5py
from scipy import stats
from colours import*
import subprocess

args = docopt(__doc__)
dir = str(args['--dir'])
transient = int(args['--t'])
fig_bool = int(args['--fig'])

## Prevent the plots ##
if fig_bool == 0:
    fig_bool = False
else:
    fig_bool = True
    
    
ThermalBoundary = []
ViscousBoundary = []
Points_in_thermal_boundary = []
Points_in_viscous_Boundary = []
run = 1
Peclet = []
Nu_integral = []
Time = []
Nu_Error = []
Re = []
Nu_top = []
Nu_bottom = []
Nu_midplane = []
u_max = []
v_max = []
w_max = []
Buoyancy = []
Dissipation = []
Energy_Balance = []
D_viscosity = []
upper_thermal_boundary = []

for idx,lines in enumerate(dir.split('/')[1].split('_')):
    if idx == 1:
        Ra = float(lines.replace('-','.'))
    if idx == 3:
        a,b,c,d,e,f,g,h = lines
        Ek = float(a+'.'+c+d+'e-'+g+h) ## This code is awful change it at some point ##
    if idx == 5:
        Pr = float(lines.replace('-','.'))
    if idx == 7:
        Nx = float(lines)  

Ny = Nx
Nz = Nx/2
Lx = Nx/Nz


onlyfiles = [f for f in os.listdir(dir) if os.path.isfile(os.path.join(dir, f))]
file = dir + '/' + onlyfiles[0]

log_file = open(file, 'r')
contents = log_file.read().split('\n')

if os.path.isdir(dir+'/img') == True:
    pass
else:
    os.system('mkdir {}/img'.format(dir))
    
    
for idx,lines in enumerate(contents):
    if idx != 0 or idx != 1:
        try:
            Time.append(float(lines.split()[0]))
            Re.append(float(lines.split()[1]))
            Nu_top.append(float(lines.split()[2]))
            Nu_bottom.append(float(lines.split()[3]))
            Nu_midplane.append(float(lines.split()[4]))
            Nu_integral.append(float(lines.split()[5]))
            u_max.append(float(lines.split()[6]))
            v_max.append(float(lines.split()[7]))
            w_max.append(float(lines.split()[8]))
            Buoyancy.append(float(lines.split()[9]))
            Dissipation.append(float(lines.split()[10]))
            Energy_Balance.append(float(lines.split()[11]))
            Nu_Error.append(float(lines.split()[12]))
            D_viscosity.append(float(lines.split()[13]))
        except:
            pass
        
Time = np.array(Time)
if fig_bool == True:

    print('The output is set to plot, creating figures in directory {}img.'.format(dir))
    Nusselt_fig = plt.figure(figsize = (12,6))
    plt.plot(Time*np.sqrt(Ra*Pr), Nu_integral, label = '$Nu_I$', color = CB91_Blue)
    plt.plot(Time*np.sqrt(Ra*Pr), Nu_bottom, label = '$Nu_b$', color = CB91_Violet)
    plt.plot(Time*np.sqrt(Ra*Pr), Nu_top, label = '$Nu_t$', color = CB91_Amber)
   # plt.plot(Time*np.sqrt(Ra*Pr), Nu_midplane, label = 'Midplane Nusselt')
    plt.xlabel('Time')
    plt.ylabel('Nusselt Number')
    plt.legend(frameon = False, ncol = 2)
    #plt.show()
    plt.savefig('{}img/Nusselt.eps'.format(dir), dpi = 500)

    Nu_error_fig = plt.figure(figsize = (10,10))
    plt.plot(Time*np.sqrt(Ra*Pr), Nu_Error)
    plt.xlabel('Time')
    plt.ylabel('Nusselt Error')
    #plt.show()
    plt.savefig('{}img/Nusselt_error.eps'.format(dir), dpi = 500)

    E_balance_fig = plt.figure(figsize = (10,10))
    plt.plot(Time*np.sqrt(Ra*Pr), Energy_Balance)
    plt.xlabel('Time')
    plt.ylabel('Energy Balance')
    #plt.show()
    plt.savefig('{}img/Energy_Balance.eps'.format(dir), dpi = 500)

    Velocity_fig = plt.figure(figsize = (10,10))
    plt.plot(Time*np.sqrt(Ra*Pr), u_max, label = 'u')
    plt.plot(Time*np.sqrt(Ra*Pr), v_max, label = 'v')
    plt.plot(Time*np.sqrt(Ra*Pr), w_max, label = 'w')
    plt.legend()
    plt.xlabel('Time')
    #plt.show()
    plt.savefig('{}img/Max_Velocity.eps'.format(dir), dpi = 500)

    Buoyancy_dissipation_fig = plt.figure(figsize = (10,10))
    plt.plot(Time*np.sqrt(Ra*Pr), Buoyancy, label = 'Buoyancy')
    plt.plot(Time*np.sqrt(Ra*Pr), Dissipation, label = 'Dissipation')
    plt.legend()
    plt.xlabel('Time')
    #plt.show()
    plt.savefig('{}img/Buoyancy_Dissipation.eps'.format(dir), dpi = 500)

    Re_fig = plt.figure(figsize = (10,10))
    plt.plot(Time*np.sqrt(Ra*Pr), Re)
    plt.xlabel('Time')
    plt.ylabel('Reynolds Number')
    #plt.show()
    plt.savefig('{}img/Reynolds.eps'.format(dir), dpi = 500)

    D_viscosity_fig = plt.figure(figsize = (10,10))
    plt.plot(Time*np.sqrt(Ra*Pr), D_viscosity)
    plt.xlabel('Time')
    plt.ylabel('Viscous Dissipation')                  
    plt.show()
    plt.savefig('{}img/D_viscosity.eps'.format(dir), dpi = 500)



    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    