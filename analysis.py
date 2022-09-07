"""Analysis file of 3d rotating rayleigh benard convection

Usage:
    analysis.py --dir=<directory> [--t=<transient> --fig=<Figure> --snap_t=<snap_t>]
    analysis.py -h | --help

Options:
    -h --help   			Display this help message
    --dir=<directory>        		Directory
    --t=<transient>   			Transient to be ignored [default: 500]
    --fig=<Figure>    			Produce Figures [default: 0]
    --snap_t=<snap_t>	                Snapshots to be ignored [default: 2]
"""

from docopt import docopt
import numpy as np
import matplotlib.pyplot as plt
import os
import h5py
from scipy import stats
from colours import*
import subprocess

# =============================================================================
# Extract the docopt arguments 
# =============================================================================

args = docopt(__doc__)
dir = str(args['--dir'])
transient = int(args['--t'])
fig_bool = int(args['--fig'])
snap_t = int(args['--snap_t'])

# =============================================================================
# Define boolean arguments for use later
# =============================================================================

if fig_bool == 0:
    fig_bool = False
else:
    fig_bool = True

# =============================================================================
# Define some empty lists which shall store data for later
# =============================================================================

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

# =============================================================================
# Extract run information
# =============================================================================

for idx,lines in enumerate(dir.split('/')[1].split('_')):
    if idx == 1:
        Ra = float(lines.replace('-','.'))
    if idx == 3:
        a,b,c,d,e,f,g,h = lines
        Ek = float(a+'.'+c+d+'e-'+g+h) 
    if idx == 5:
        Pr = float(lines.replace('-','.'))
    if idx == 7:
        Nx = float(lines)  

Ny = Nx
Nz = Nx/2
Lx = Nx/Nz

# =============================================================================
# Find the name of the log file ans store it in fileName
# =============================================================================

onlyfiles = [f for f in os.listdir(dir) if os.path.isfile(os.path.join(dir, f))]
fileName = dir + '/' + onlyfiles[0]

# =============================================================================
# Extract the data from the log file
# =============================================================================

logFile = open(fileName, 'r')
logFileContents = logFile.read().split('\n')

# =============================================================================
# Check if an image directory already exsists, if not create one
# =============================================================================

if os.path.isdir(dir+'/img') == True:
    pass
else:
    os.system('mkdir {}/img'.format(dir))
    
# =============================================================================
# Fill up the arrays with data from the log file
# =============================================================================

for idx,lines in enumerate(logFileContents):
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

# =============================================================================
# Extract the energy balance from the log file
# =============================================================================

E_B = np.array(Energy_Balance)
E_B[E_B == np.inf] = 0
Energy_Balance = E_B[np.logical_not(np.isnan(E_B))]
Balance = np.average(Energy_Balance[-transient:])*100
logFile.close()

# =============================================================================
# Deal with the analysis tasks to calculate boundary layers
# =============================================================================

analysis_file = os.listdir(dir)

# =============================================================================
# Fine the file names and store them as a variable 
# =============================================================================

for idx,lines in enumerate(os.listdir(dir)):
    if lines.find('analysis') != -1:
        analysis_file_name = os.listdir(dir)[idx]

for idx,lines in enumerate(os.listdir(dir)):
    if lines.find('snapshot') != -1:
        snapshot_file_name = os.listdir(dir)[idx]

# =============================================================================
# Load the data from the analysis file 
# =============================================================================

with h5py.File('{}{}/analysis.h5'.format(dir,analysis_file_name), mode = 'r') as file:
    U_h = np.copy(file['tasks']['U_H_prof'])
    z = np.copy(file['tasks']['z'])[-1,0,0,:]
    conduction = np.copy(file['tasks']["conduction_prof"])
    advection = np.copy(file['tasks']["advection_prof"])
    Pe = np.copy(file['tasks']['Pe'])
    T_prof = np.copy(file['tasks']["T_prof"])    
    thermal_dissipation = np.copy(file['tasks']["thermal_dissipation"])[-transient:,0,0,0] 
    
# =============================================================================
# Load the data from the snapshot file 
# =============================================================================

with h5py.File('{}{}/snapshots.h5'.format(dir,snapshot_file_name), mode = 'r') as file:
   full_T = np.copy(file['tasks']['T'])
   T = full_T[-snap_t:,:,:,:]


def derivative(data,z):
    """
    

    Parameters
    ----------
    data : TYPE
        DESCRIPTION.
    z : TYPE
        DESCRIPTION.

    Returns
    -------
    The derivative of the data with respect to z.

    """
    return np.gradient(data,z)

def determine_root(derivative,z):
    """
    

    Parameters
    ----------
    derivative : TYPE
        DESCRIPTION.
    z : TYPE
        DESCRIPTION.

    Returns
    -------
    upper : The upper boundary thickness
    lower : The lower boundary thickness
    avg : The average boundary thickness
    avgPoints : Average number of points in the boundary

    """

    upper,lower,avg = 0,0,0
    zero_crossings = np.where(np.diff(np.sign(derivative)))[0]

    # =========================================================================
    # Determine the lower crossing
    # =========================================================================
    
    x,y = [z[zero_crossings[0]],z[zero_crossings[0]+1]],[derivative[zero_crossings[0]],derivative[zero_crossings[0]+1]]
    m,b = np.polyfit(y,x,1)
    lower = b

    # =========================================================================
    # Determine the upper crossing 
    # =========================================================================

    x,y = [z[zero_crossings[-1]],z[zero_crossings[-1]+1]],[derivative[zero_crossings[-1]],derivative[zero_crossings[-1]+1]]
    m,b = np.polyfit(y,x,1)
    upper = b
    avg = (lower + (1-upper))/2

    # =========================================================================
    # Calculate avg points in boundary 
    # =========================================================================

    avgPoints = zero_crossings[0]+1
    return upper, lower, avg, avgPoints

# =============================================================================
# Construct a 3d array the same size as T in order to use numpy to subtract 
# T_prof, i.e. avoid for loops. - Not sure if this is still needed, will check
# on arc.
# =============================================================================

avg_T_prof = np.average(np.array(T_prof[:,0,0,:]), axis=0) 
T_rms_prof = np.zeros(len(z))

for idxt,time_step in enumerate(T):
    if idxt >= 0:
        for idx,lines in enumerate(np.rot90(time_step, k = 3, axes = (0,2))):
            T_rms_prof[idx] += np.average(np.sqrt((lines - avg_T_prof[idx]*np.ones(np.shape(lines)))**2))

T_rms_derivative = derivative(T_rms_prof, z)
upper_trms,lower_trms,avg_trms,points_trms = determine_root(T_rms_derivative, z)

# =============================================================================
# Calculating the thermal boundary layer, should probably re-write at some 
# point. Note great code.
# =============================================================================

for idx,lines in enumerate(advection):
    if idx > np.shape(conduction)[0] - transient:
        try:
            profile = (advection[idx,0,0,:] + conduction[idx,0,0,:])
            idx1 = np.where(profile[:-1] * profile[1:] < 0 )[0] + 1
            Points_in_thermal_boundary.append(idx1[0] + 1)
            upper_z = [z[idx1[-1]],z[idx1[-1]+1]]
            upper_profile = [float(profile[idx1[-1]]), float(profile[idx1[0]+1])]
            zz = [z[idx1[0]-1],z[idx1[0]]]
            pprofile = [float(profile[idx1[0] - 1]), float(profile[idx1[0]])]
            upper_m, upper_b = np.polyfit(upper_profile, upper_z, 1)
            upper_thermal_boundary.append(0.5-upper_b)
            m, b = np.polyfit(pprofile, zz, 1)
            ThermalBoundary.append(abs(0.5+b))
        except:
            pass
        
# =============================================================================
# Calculate the average thermal boundary layer thickness using Rob Longs method
# =============================================================================

avg_temp_profile = np.average(np.array(advection[-transient:,0,0,:] + conduction[-transient:,0,0,:]), axis=0)

# =============================================================================
# Avg the horizontal velocity profile over time 
# =============================================================================
avg_u_h = np.average(np.array(U_h[-transient:,0,0,:]), axis=0)
upper_viscous_boundary, lower_viscous_boundary, avg_viscous_boundary, avg_points = determine_root(derivative(avg_u_h,z),z)

Dt = 0.5 + lower_trms
dv = avg_viscous_boundary
points_in_dt = points_trms
points_in_dv = avg_points

# =============================================================================
# Print everything to the terminal 
# =============================================================================

widthOfTerminal = int(subprocess.check_output("tput cols", shell=True))
print('\n')
print('-'*widthOfTerminal)
print('|     -    | Run time  |   Convective time | File length | Snapshots |     Pe    |    Re    |  Int Nu  |   Err Nu   |  Viscous/Points  | Thermal/Points | Energy Balance |')
print('-'*widthOfTerminal)
print('| Avgerage |   {:.2f}    |       {}        |    {}    |     {}    |  {:.3f}  |  {:.3f}  |   {:.3f}  |     {:.2f}   |   {:.3f}/{:.1f}   |   {:.3f}/{:.1f} |      {:.2f}%     |'.format(Time[-1], int(Time[-1]*np.sqrt(Ra*Pr)), len(Nu_integral), len(full_T), np.average(Pe[-transient:]), np.average(Re[-transient:]), np.average(Nu_integral[-transient:]), np.average(Nu_Error[-transient:]), dv, avg_points, Dt, points_trms, Balance))
print('-'*widthOfTerminal)
print('|   std    |     -     |        -          |      -      |     -     |   {:.3f}   |   {:.3f}  |   {:.3f}  |     {:.2f}   |     -     |     -   |        -       |'.format(np.std(Pe[-transient:]), np.std(Re[-transient:]), np.std(Nu_integral[-transient:]), np.std(Nu_Error[-transient:])     ))
print('-'*widthOfTerminal)
print('\n')
print('Buoyancy = {:.4e}, Dissipation = {:.4e}, Energy Balance = {:.3f}%.'.format(np.average(Buoyancy[-transient:]), np.average(Dissipation[-transient:]), Balance))

# =============================================================================
# print('Top Nusselt: {:.3f}, std = {:.3f}.'.format(np.average(Nu_top[-transient:]),np.std(Nu_top[-transient:])))
# print('Bottom Nusselt: {:.3f}, std = {:.3f}.'.format(np.average(Nu_bottom[-transient:]),np.std(Nu_bottom[-transient:])))
# print('Midplane Nusselt: {:.3f}, std = {:.3f}.'.format(np.average(Nu_midplane[-transient:]),np.std(Nu_midplane[-transient:])))
# print('Max Velocities: u = {:.3f}, v = {:.3f}, w = {:.3f}.'.format(np.average(u_max), np.average(v_max), np.average(w_max)))
# print('Using (Long,2020) upper: {:.4f}, lower: {:.4f}, avg: {:.4f}, points in boundary: {}.'.format(np.average(ThermalBoundary),np.average(upper_thermal_boundary),(np.average(ThermalBoundary)+np.average(upper_thermal_boundary))/2, points_in_dt))
# print('The average viscous boundary layer thickness is {:.4f}, and there are {} points in the boundary.'.format(avg_viscous_boundary, avg_points))
# print('Using the T rms method, upper: {:.5f}, lower: {:.5f}, avg: {:.5f} and number of points {}.'.format(0.5-upper_trms,0.5+lower_trms,avg_trms,points_trms))
# =============================================================================

# =============================================================================
# Convert time array to numpy for manipulation 
# =============================================================================

Time = np.array(Time) 

if fig_bool == True:

    print('The output is set to plot, creating figures in directory {}img.'.format(dir))
    Nusselt_fig = plt.figure(figsize = (12,6))
    plt.plot(Time*np.sqrt(Ra*Pr), Nu_integral, label = '$Nu_I$', color = CB91_Blue)
    plt.plot(Time*np.sqrt(Ra*Pr), Nu_bottom, label = '$Nu_b$', color = CB91_Violet)
    plt.plot(Time*np.sqrt(Ra*Pr), Nu_top, label = '$Nu_t$', color = CB91_Amber)
    plt.xlabel('Time')
    plt.ylabel('Nusselt Number')
    plt.legend(frameon = False, ncol = 2)
    plt.savefig('{}img/Nusselt.eps'.format(dir), dpi = 500)

    Nu_error_fig = plt.figure(figsize = (10,10))
    plt.plot(Time*np.sqrt(Ra*Pr), Nu_Error)
    plt.xlabel('Time')
    plt.ylabel('Nusselt Error')
    plt.savefig('{}img/Nusselt_error.eps'.format(dir), dpi = 500)

    E_balance_fig = plt.figure(figsize = (10,10))
    plt.plot(Time*np.sqrt(Ra*Pr), Energy_Balance)
    plt.xlabel('Time')
    plt.ylabel('Energy Balance')
    plt.savefig('{}img/Energy_Balance.eps'.format(dir), dpi = 500)

    Velocity_fig = plt.figure(figsize = (10,10))
    plt.plot(Time*np.sqrt(Ra*Pr), u_max, label = 'u')
    plt.plot(Time*np.sqrt(Ra*Pr), v_max, label = 'v')
    plt.plot(Time*np.sqrt(Ra*Pr), w_max, label = 'w')
    plt.legend()
    plt.xlabel('Time')
    plt.savefig('{}img/Max_Velocity.eps'.format(dir), dpi = 500)

    Buoyancy_dissipation_fig = plt.figure(figsize = (10,10))
    plt.plot(Time*np.sqrt(Ra*Pr), Buoyancy, label = 'Buoyancy')
    plt.plot(Time*np.sqrt(Ra*Pr), Dissipation, label = 'Dissipation')
    plt.legend()
    plt.xlabel('Time')
    plt.savefig('{}img/Buoyancy_Dissipation.eps'.format(dir), dpi = 500)

    Re_fig = plt.figure(figsize = (10,10))
    plt.plot(Time*np.sqrt(Ra*Pr), Re)
    plt.xlabel('Time')
    plt.ylabel('Reynolds Number')
    plt.savefig('{}img/Reynolds.eps'.format(dir), dpi = 500)

    D_viscosity_fig = plt.figure(figsize = (10,10))
    plt.plot(Time*np.sqrt(Ra*Pr), D_viscosity)
    plt.xlabel('Time')
    plt.ylabel('Viscous Dissipation')                  
    plt.savefig('{}img/D_viscosity.eps'.format(dir), dpi = 500)


    avg_u_h_viscosity_fig = plt.figure(figsize = (10,10))
    plt.plot(z, avg_u_h)
    plt.axvline(lower_viscous_boundary, color = 'black', lw = 2)
    plt.axvline(upper_viscous_boundary, color = 'black', lw = 2)
    plt.xlabel('z')
    plt.ylabel('$U_h$')
    plt.savefig('{}img/avg_u_h.eps'.format(dir), dpi = 500)


    temp_boundary_fig = plt.figure(figsize = (10,10))
    plt.plot(z, avg_temp_profile)
    plt.xlabel('z')
    plt.ylabel('Temperature Profile')
    plt.savefig('{}img/temp_boundary.eps'.format(dir), dpi = 500)


    temp_rms_fig = plt.figure(figsize = (10,10))
    plt.plot(z, T_rms_prof)
    plt.axvline(lower_trms, color = 'black', lw = 2)
    plt.axvline(upper_trms, color = 'black', lw = 2)
    plt.plot(z,np.zeros(len(z)))
    plt.xlabel('z')
    plt.ylabel('$T_{rms}$')
    plt.savefig('{}img/trms_boundary.eps'.format(dir), dpi = 500)