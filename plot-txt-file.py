"""Analysis file of 3d rotating rayleigh benard convection

Usage:
    plot-txt-file.py --dir=<directory> [--t=<transient>]
    plot-txt-file.py -h | --help

Options:
    -h --help   		  Display this help message
    --dir=<directory>             Directory
    --t=<transient>               Transient to be ignored [default: 2000]
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
lineWidth = 2


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


onlyFiles = [f for f in os.listdir(dir) if os.path.isfile(os.path.join(dir, f))]

for lines in onlyFiles:

    if lines.find('.txt') != 1:
        logFile = dir + '/' + lines

log_file = open(logFile, 'r')
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
iterations = len(Time)*10
lastSnapShot = iterations - (iterations % 5000) + 5000
snapShots = [Time[int(snaps/10)]*np.sqrt(Ra*Pr) for snaps in range(lastSnapShot - 25000, lastSnapShot, 5000)]

for i in range(lastSnapShot - 25000, lastSnapShot, 5000):
    print(i)

print('The output is set to plot, creating figures in directory {}img.'.format(dir))
Nusselt_fig = plt.figure(figsize = (12,6))
plt.plot(Time*np.sqrt(Ra*Pr), Nu_integral, label = '$Nu_I$', color = CB91_Blue, lw = lineWidth)
plt.plot(Time*np.sqrt(Ra*Pr), Nu_bottom, label = '$Nu_b$', color = CB91_Violet, lw = lineWidth)
plt.plot(Time*np.sqrt(Ra*Pr), Nu_top, label = '$Nu_t$', color = CB91_Amber, lw = lineWidth)

for snaps in snapShots:
    plt.plot(np.ones(100)*snaps, np.linspace(0, max(Nu_integral), 100),'k--', lw = 1)

plt.xlabel('Time')
plt.ylabel('Nusselt Number')
plt.title('Nusselt Number Comparison')
plt.legend(frameon = False, ncol = 2)
plt.savefig('{}img/Nusselt.eps'.format(dir), dpi = 500)

Velocity_fig = plt.figure(figsize = (12,6))
plt.plot(Time*np.sqrt(Ra*Pr), u_max, label = 'u', color = CB91_Blue, lw = lineWidth)
plt.plot(Time*np.sqrt(Ra*Pr), v_max, label = 'v', color = CB91_Violet, lw = lineWidth)
plt.plot(Time*np.sqrt(Ra*Pr), w_max, label = 'w', color = CB91_Amber, lw = lineWidth)

for snaps in snapShots:
    plt.plot(np.ones(100)*snaps, np.linspace(0, max(u_max), 100),'k--', lw = 1)

plt.legend(frameon = False, ncol = 2)
plt.xlabel('Time')
plt.savefig('{}img/Max_Velocity.eps'.format(dir), dpi = 500)

Buoyancy_dissipation_fig = plt.figure(figsize = (12,6))
plt.plot(Time*np.sqrt(Ra*Pr), Buoyancy, label = 'Buoyancy', lw = lineWidth, color = CB91_Blue)
plt.plot(Time*np.sqrt(Ra*Pr), Dissipation, label = 'Dissipation', lw = lineWidth, color = CB91_Amber)

for snaps in snapShots:
    plt.plot(np.ones(100)*snaps, np.linspace(0, max(Buoyancy), 100),'k--', lw = 1)

plt.legend(frameon = False, ncol = 2)
plt.xlabel('Time')
plt.savefig('{}img/Buoyancy_Dissipation.eps'.format(dir), dpi = 500)

Re_fig = plt.figure(figsize = (12,6))

for snaps in snapShots:
    plt.plot(np.ones(100)*snaps, np.linspace(0, max(Nu_integral), 100),'k--', lw = 1)

plt.plot(Time*np.sqrt(Ra*Pr), Re, lw = lineWidth, color = CB91_Blue)
plt.xlabel('Time')
plt.ylabel('Reynolds Number')
plt.savefig('{}img/Reynolds.eps'.format(dir), dpi = 500)

D_viscosity_fig = plt.figure(figsize = (12,6))
plt.plot(Time*np.sqrt(Ra*Pr), D_viscosity, lw = lineWidth, color = CB91_Blue)

for snaps in snapShots:
    plt.plot(np.ones(100)*snaps, np.linspace(0, max(D_viscosity), 100),'k--', lw = 1)

plt.xlabel('Time')
plt.ylabel('Viscous Dissipation')                  
plt.show()
plt.savefig('{}img/D_viscosity.eps'.format(dir), dpi = 500)



    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
