"""Analysis of forces file of 3d rotating rayleigh benard convection

Usage:
    forces_spectrums.py --dir=<directory> --snap_t=<snapshot_transient> --mask=<mask> [--t=<transient> --fig=<Figure>]
    forces_spectrums.py -h | --help

Options:
    -h --help                           Display this help message
    --dir=<directory>                   Directory
    --snap_t=<transient>                Snapshot transient
    --mask=<mask>                       Number of viscous boundaries to ignore 
    --t=<transient>                     Transient to be ignored [default: 2000]
    --fig=<Figure>                      Produce Figures [default: True]
"""

from docopt import docopt
import numpy as np
import matplotlib.pyplot as plt
import os
import h5py
from scipy import stats
from dedalus import public as de
from dedalus.core.operators import Integrate 
from colours import *

args = docopt(__doc__)
dir = str(args['--dir'])
transient = int(args['--t'])
fig_bool = bool(args['--fig'])
transient = int(args['--t'])
snap_t = int(args['--snap_t'])
mask = int(args['--mask'])

plt.rcParams['font.family'] = 'Serif'
plt.rcParams['font.size'] = 18
plt.rcParams['axes.linewidth'] = 1

BuoyancyColour = CB91_Green
ViscosityColour = CB91_Pink
InertiaColour = 'maroon'
CoriolisColour = CB91_Violet
PressureColour = 'darkorange'
ACColour = CB91_Blue
spectrumlw = 2


for idx,lines in enumerate(os.listdir(dir)):
    if lines.find('analysis') != -1:
        forces_file_name = os.listdir(dir)[idx]

for idx,lines in enumerate(dir.split('/')[1].split('_')):
    if idx == 1:
        Ra = float(lines.replace('-','.'))
    if idx == 3:
        a,b,c,d,e,f,g,h = lines
        Ek = float(a+'.'+c+d+'e-'+g+h) ## This code is awful change it at some point ##
    if idx == 5:
        Pr = float(lines.replace('-','.'))
    if idx == 7:
        N = int(lines)

with h5py.File('{}{}/analysis.h5'.format(dir,forces_file_name), mode = 'r') as file:

    z = np.copy(file['tasks']["z"])[0,0,0,:]
    U_h = np.copy(file['tasks']['U_H_prof'])

if os.path.isdir(dir+'/img') == True:
    pass
else:
    os.system('mkdir {}/img'.format(dir))


for idx,lines in enumerate(os.listdir(dir)):
    if lines.find('snapshot') != -1:
        snapshot_file_name = os.listdir(dir)[idx]


with h5py.File('{}{}/snapshots.h5'.format(dir,snapshot_file_name), mode = 'r') as file:

    ## Velocity data ##
    kinetic_spectrum = np.copy(file['tasks']["kinetic_spectrum"])[-snap_t:,:,:,:]
    u = np.copy(file['tasks']["u"])[-snap_t:,:,:,:]
    v = np.copy(file['tasks']["v"])[-snap_t:,:,:,:]
    w = np.copy(file['tasks']["w"])[-snap_t:,:,:,:]
    T = np.copy(file['tasks']["T"])[-snap_t:,:,:,:]

    ## X equation
    x_pressure = np.copy(file['tasks']["x_pressure"])[-snap_t:,:,:,:]
    x_diffusion = np.copy(file['tasks']["x_diffusion"])[-snap_t:,:,:,:]
    x_coriolis = np.copy(file['tasks']["x_coriolis"])[-snap_t:,:,:,:]
    x_inertia = np.copy(file['tasks']["x_inertia"])[-snap_t:,:,:,:]

    ## Y equation
    y_pressure = np.copy(file['tasks']["y_pressure"])[-snap_t:,:,:,:]
    y_diffusion = np.copy(file['tasks']["y_diffusion"])[-snap_t:,:,:,:]
    y_coriolis = np.copy(file['tasks']["y_coriolis"])[-snap_t:,:,:,:]
    y_inertia =	np.copy(file['tasks']["y_inertia"])[-snap_t:,:,:,:]

    ## Z equation
    z_pressure = np.copy(file['tasks']["z_pressure"])[-snap_t:,:,:,:]
    z_diffusion = np.copy(file['tasks']["z_diffusion"])[-snap_t:,:,:,:]
    z_inertia = np.copy(file['tasks']["z_inertia"])[-snap_t:,:,:,:]
    z_buoyancy = np.copy(file['tasks']["z_bouyancy"])[-snap_t:,:,:,:]


####################################################################################################
## Some code to try and take the spectrum of the real fields, i.e. get the spectrum of the forces ##
####################################################################################################

Nx = N 
Ny = N
Lx = Ly = 1
Lz = 1
Nz = int(N/Lx)

x_basis = de.Fourier('x', Nx, interval = (0,Lx), dealias=3/2)
y_basis = de.Fourier('y', Ny, interval = (0,Ly), dealias=3/2)
z_basis = de.Chebyshev('z', Nz, interval = (-Lz/2,Lz/2), dealias =3/2)
domain = de.Domain([x_basis, y_basis, z_basis], grid_dtype=np.float64)


## Create fields in dedalus so we can use dedalus tools for computations, retain spectral accuracy e.c.t ##
Viscosity = domain.new_field(name='viscosity')
Coriolis = domain.new_field(name='coriolis')
ACoriolis = domain.new_field(name='Acoriolis')
Inertia = domain.new_field(name='inertia')
Buoyancy = domain.new_field(name='buoyancy')
Pressure = domain.new_field(name='pressure')
Kinetic = domain.new_field(name='kinetic')

print(np.shape(Viscosity['c']))

## Velocity fields ##
U = domain.new_field(name='u')
V = domain.new_field(name='v')
W = domain.new_field(name='w')
Nusselt = domain.new_field(name='Nusselt')

def ComputeRMS(Fx, Fy, Fz):

    ##################################################################################
    ## Like the equation used in (Andres, 2021) this computes the RMS of the force. ##
    ##################################################################################

    return np.sqrt(Fx**2 + Fy**2 + Fz**2)


def Compute(Force):

    ##################################################################################
    ## Function which computes all the usefull information that we want to extract, ##
    ## takes the argument 'Force' which is the force we are using.                  ##
    ##################################################################################

    #Force.meta['z']['constant'] = True
    
    ForceSpectrum = (Force['c'].imag**2 +  Force['c'].real**2)**0.5
    z_avg = np.sum(ForceSpectrum, axis = 0) # average over z 
    return np.sum(z_avg, axis=1) # average over x or y, change the axis to pick which one

def HorizontalAverage(ForceRMS):

    #####################################################################################################################################
    ## Takes in the rms of the force and returns a 1d array which contains the plane averaged profile, i.e. F(z) as in (Guzman, 2021). ##
    #####################################################################################################################################
    profile = []

    RotatedForce = np.rot90(ForceRMS, k=1, axes=(2,0))
 
    for lines in RotatedForce:
        profile.append(np.sum(lines))

    return np.array(profile)

def derivative(data,z):
    derivative = np.gradient(data,z)
    return derivative

def determine_root(derivative,z):

    upper,lower,avg = 0,0,0
    zero_crossings = np.where(np.diff(np.sign(derivative)))[0]

    ## Determine the lower crossing ##
    x,y = [z[zero_crossings[0]],z[zero_crossings[0]+1]],[derivative[zero_crossings[0]],derivative[zero_crossings[0]+1]]
    m,b = np.polyfit(y,x,1)
    lower = b

    ## Determine the upper crossing ##
    x,y = [z[zero_crossings[-1]],z[zero_crossings[-1]+1]],[derivative[zero_crossings[-1]],derivative[zero_crossings[-1]+1]]
    m,b = np.polyfit(y,x,1)
    upper = b
    avg = (lower + (1-upper))/2

    ## Calculate avg points in boundary ##
    avg_points = zero_crossings[0]*mask

    return upper, lower, avg, avg_points


def RemoveBoundaries(Mask, Force):

    RotatedForce = np.rot90(Force['g'], k=1, axes = (2,0))
    MaskedForce = RotatedForce*Mask
    Force['g'] = np.rot90(MaskedForce,	k=-1, axes = (2,0))

## Avg the horizontal velocity profile over time ##
avg_u_h = np.average(np.array(U_h[-transient:,0,0,:]), axis=0) ## Create an array which contains only the last N points, then avg over this ##
upper_viscous_boundary, lower_viscous_boundary, avg_viscous_boundary, avg_points = determine_root(derivative(avg_u_h,z),z)


## Construct a matrix which is 1 everywhere other than in the viscous boundary where it is 0. This can be used to remove the viscous boundaries 
## from my force plots

Mask = np.rot90(np.ones(np.shape(x_diffusion[-1]),dtype=np.int32), k=1 , axes=(2,0)) ## Rotate that boi


for idx,slices in enumerate(Mask):
    
    if idx < avg_points or idx >= len(z) - avg_points:

        Mask[idx] = np.zeros(np.shape(slices))

## Arrays containing the time series ##
ViscosityTimeSeries = []
CoriolisTimeSeries = []
BuoyancyTimeSeries = []
InertiaTimeSeries = []
PressureTimeSeries = []
ACoriolisTimeSeries = []
KineticTimeSeries = []
UTimeSeries = []
VTimeSeries = []
WTimeSeries = []


## Arrays containing the time series ##
HAvgViscosityTimeSeries = []
HAvgCoriolisTimeSeries = []
HAvgBuoyancyTimeSeries = []
HAvgInertiaTimeSeries = []
HAvgPressureTimeSeries = []
HAvgACoriolisTimeSeries = []
BlankMatrix = np.zeros(np.shape(z_diffusion[-1]), dtype=np.int32)

for idx in range(1,snap_t+1):

    ## Load data
    Viscosity['g'] = ComputeRMS(x_diffusion[-idx], y_diffusion[-idx], z_diffusion[-idx])
    Coriolis['g'] = ComputeRMS(x_coriolis[-idx], y_coriolis[-idx], BlankMatrix)
    Inertia['g'] = ComputeRMS(x_inertia[-idx], y_inertia[-idx], z_inertia[-idx])
    Buoyancy['g'] = ComputeRMS(BlankMatrix, BlankMatrix, z_buoyancy[-idx])
    Pressure['g'] = ComputeRMS(x_pressure[-idx], y_pressure[-idx], z_pressure[-idx])
    Kinetic['c'] = kinetic_spectrum[-idx]
    ACoriolis['g'] = ComputeRMS(x_pressure[-idx] + x_coriolis[-idx], y_pressure[-idx] + y_coriolis[-idx], z_pressure[-idx] - BlankMatrix)
    U['g'] = u[-idx]
    V['g'] = v[-idx]
    W['g'] = w[-idx]
    
    ## Remove Boundaries
    RemoveBoundaries(Mask, Viscosity)
    RemoveBoundaries(Mask, Pressure)
    RemoveBoundaries(Mask, Coriolis)
    RemoveBoundaries(Mask, Inertia)
    RemoveBoundaries(Mask, ACoriolis)
    RemoveBoundaries(Mask, Buoyancy)
    RemoveBoundaries(Mask, Kinetic)
    RemoveBoundaries(Mask, U)
    RemoveBoundaries(Mask, V)
    RemoveBoundaries(Mask, W)


    
    ## Compute the spectrums of each force ##
    ViscosityTimeSeries.append(Compute(Viscosity))
    CoriolisTimeSeries.append(Compute(Coriolis))
    InertiaTimeSeries.append(Compute(Inertia))
    BuoyancyTimeSeries.append(Compute(Buoyancy))
    PressureTimeSeries.append(Compute(Pressure))
    KineticTimeSeries.append(Compute(Kinetic))
    UTimeSeries.append(Compute(U))
    VTimeSeries.append(Compute(V))
    WTimeSeries.append(Compute(W))
    ACoriolisTimeSeries.append(Compute(ACoriolis))

    ## Compute the horizontal average ##
    HAvgViscosityTimeSeries.append(HorizontalAverage(Viscosity['g']))
    HAvgCoriolisTimeSeries.append(HorizontalAverage(Coriolis['g']))
    HAvgInertiaTimeSeries.append(HorizontalAverage(Inertia['g']))
    HAvgBuoyancyTimeSeries.append(HorizontalAverage(Buoyancy['g']))
    HAvgPressureTimeSeries.append(HorizontalAverage(Pressure['g']))
    HAvgACoriolisTimeSeries.append(abs(HorizontalAverage(ACoriolis['g'])))

    print('Coriolis: {:.2e}, Pressure: {:.2e}, Viscosity: {:.2e}, Inertia: {:.2e}, Buoyancy: {:.2e}, Ageostrophic: {:.2e}'.format(np.sum(Coriolis['g']), np.sum(Pressure['g']),np.sum(Viscosity['g']), np.sum(Inertia['g']), np.sum(Buoyancy['g']), np.sum(ACoriolis['g'])))

## Time avg the spectrums ##
ViscositySpectrum = np.average(np.array(ViscosityTimeSeries), axis=0)
InertiaSpectrum = np.average(np.array(InertiaTimeSeries), axis=0)
BuoyancySpectrum = np.average(np.array(BuoyancyTimeSeries), axis=0)
CoriolisSpectrum = np.average(np.array(CoriolisTimeSeries), axis=0)
PressureSpectrum = np.average(np.array(PressureTimeSeries), axis=0)
ACoriolisSpectrum = np.average(np.array(ACoriolisTimeSeries), axis=0)
KineticSpectrum = np.average(np.array(KineticTimeSeries), axis=0)
USpectrum = np.average(np.array(UTimeSeries), axis=0)
VSpectrum = np.average(np.array(VTimeSeries), axis=0)
WSpectrum = np.average(np.array(WTimeSeries), axis=0)

## Time avg the profiles ##
ViscosityProfile = np.average(np.array(HAvgViscosityTimeSeries), axis=0)
InertiaProfile = np.average(np.array(HAvgInertiaTimeSeries), axis=0)
BuoyancyProfile = np.average(np.array(HAvgBuoyancyTimeSeries), axis=0)
CoriolisProfile = np.average(np.array(HAvgCoriolisTimeSeries), axis=0)
PressureProfile = np.average(np.array(HAvgPressureTimeSeries), axis=0)
ACoriolisProfile = np.average(np.array(HAvgACoriolisTimeSeries), axis=0)


fig = plt.figure(figsize=(10,10))
plt.plot(range(len(ViscositySpectrum)), ViscositySpectrum, label = '$F_v$', color = ViscosityColour, lw=spectrumlw)
plt.plot(range(len(ViscositySpectrum)), CoriolisSpectrum, label = '$F_C$', color = CoriolisColour, lw=spectrumlw)
plt.plot(range(len(InertiaSpectrum)), InertiaSpectrum, label = '$F_I$', color = InertiaColour, lw=spectrumlw)
plt.plot(range(len(BuoyancySpectrum)), BuoyancySpectrum, label = '$F_B$', color = BuoyancyColour, lw=spectrumlw)
plt.plot(range(len(BuoyancySpectrum)), PressureSpectrum, label = '$F_P$', color = PressureColour, lw=spectrumlw)
plt.plot(range(len(ViscositySpectrum)), ACoriolisSpectrum, label = '$F_{AC}$', color = ACColour, lw=spectrumlw)
plt.xscale("log")
plt.yscale("log")
plt.xlabel('$K_z$')
plt.ylabel('Magnitude')
plt.xlim(1,64)
legend_properties = {'weight':'bold'}
plt.legend(ncol=2, fontsize=14,prop=legend_properties,frameon=False)
plt.savefig('{}/img/ForceSpectrum.eps'.format(dir), dpi=500)
plt.show()


print('The max wavenumber for k.e: {}, u: {}, v: {}, w: {}.'.format(np.argmax(KineticSpectrum[:64]), np.argmax(USpectrum[:64]), np.argmax(VSpectrum[:64]), np.argmax(WSpectrum[:64])))
fig = plt.figure(figsize=(12,6))
plt.plot(range(len(ViscositySpectrum)), KineticSpectrum, lw = 2, color = ViscosityColour, label = 'K.e')
plt.plot(range(len(ViscositySpectrum)), USpectrum, lw = 2, color = CoriolisColour, label = 'u')
plt.plot(range(len(ViscositySpectrum)), VSpectrum, lw = 2, color = PressureColour, label = 'v')
plt.plot(range(len(ViscositySpectrum)), WSpectrum, lw = 2, color = ACColour, label = 'w')
plt.xscale("log")
plt.yscale("log")
plt.xlabel('$K_z$')
plt.ylabel('Magnitude')
plt.xlim(1,64)
legend_properties = {'weight':'bold'}
plt.legend(ncol=2, fontsize=14,prop=legend_properties,frameon=False)
plt.savefig('{}/img/KineticSpectrum.eps'.format(dir), dpi=500)
plt.show()


## Plot the horizontal profiles ##
fig = plt.figure(figsize=(10,10))
plt.plot(ViscosityProfile, z, color = ViscosityColour, label = '$F_V$', lw = 2)
plt.plot(CoriolisProfile, z, color = CoriolisColour, label = '$F_C$', lw = 2)
plt.plot(InertiaProfile, z, color = InertiaColour, label = '$F_I$', lw = 2)
plt.plot(ACoriolisProfile, z, color = ACColour, label = '$F_{AC}$', lw = 2)
plt.plot(BuoyancyProfile, z, color = BuoyancyColour, label = '$F_B$', lw = 2)
plt.plot(PressureProfile, z, color = PressureColour, label = '$F_P$', lw = 2)
legend_properties = {'weight':'bold'}
plt.legend(ncol=2, fontsize=14,prop=legend_properties,frameon=False)
plt.xlabel('Magnitude')
plt.ylabel('z')
plt.savefig('{}/img/ForceProfiles.eps'.format(dir), dpi=500)
plt.show()





## Plot mid height values like (Guzman, 2021) ##

mid_point = int(len(z)/2)

fig = plt.figure(figsize=(12,6))
plt.xscale("log")
plt.yscale("log")
plt.scatter([Ra], ViscosityProfile[mid_point], label = 'Viscosity')
plt.scatter([Ra], CoriolisProfile[mid_point],  label = 'Coriolis')
plt.scatter([Ra], InertiaProfile[mid_point], label = 'Inertia')
plt.scatter([Ra], BuoyancyProfile[mid_point], label = 'Buoyancy')
plt.scatter([Ra], PressureProfile[mid_point], label = 'Pressure')
plt.legend()
plt.show()

