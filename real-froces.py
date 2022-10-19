"""Analysis of forces file of 3d rotating rayleigh benard convection

Usage:
    forces_spectrums.py [--dir=<directory>] [--snap_t=<snapshot_transient>] [--mask=<mask>] [--t=<transient> --fig=<Figure>]
    forces_spectrums.py -h | --help

Options:
    -h --help                           Display this help message
    --dir=<directory>                   Directory [default: results/Ra_1-40e+08_Ek_1-00e-05_Pr_7-0_N_128_q_1-2_k_48_enhanced/]
    --snap_t=<transient>                Snapshot transient [default: 1]
    --mask=<mask>                       Number of viscous boundaries to ignore [default: 0] 
    --t=<transient>                     Transient to be ignored [default: 2000]
    --fig=<Figure>                      Produce Figures [default: False]
"""

from docopt import docopt
import numpy as np
import matplotlib.pyplot as plt
import os
import h5py
from scipy import stats
from scipy.fft import fft
from dedalus import public as de
from dedalus.core.operators import Integrate 
from colours import *

%matplotlib inline
# =============================================================================
# Extract the docopt arguments 
# =============================================================================

args = docopt(__doc__)
dir = str(args['--dir'])
transient = int(args['--t'])
fig_bool = bool(args['--fig'])
transient = int(args['--t'])
snap_t = int(args['--snap_t'])
mask = int(args['--mask'])

# =============================================================================
# Set some plotting parameters for later 
# =============================================================================

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

# =============================================================================
# Create the string which leads to the right directory
# =============================================================================

for idx,lines in enumerate(os.listdir(dir)):
    if lines.find('analysis') != -1:
        forces_file_name = os.listdir(dir)[idx]

for idx,lines in enumerate(dir.split('/')[1].split('_')):
    if idx == 1:
        Ra = float(lines.replace('-','.'))
    if idx == 3:
        a,b,c,d,e,f,g,h = lines
        Ek = float(a+'.'+c+d+'e-'+g+h)
    if idx == 5:
        Pr = float(lines.replace('-','.'))
    if idx == 7:
        N = int(lines)
        
# =============================================================================
# Extract some profiles from the analysis file, these are used later
# =============================================================================

with h5py.File('{}{}/analysis.h5'.format(dir,forces_file_name), mode = 'r') as file:

    z = np.copy(file['tasks']["z"])[0,0,0,:]
    U_h = np.copy(file['tasks']['U_H_prof'])
    
# =============================================================================
# Check to see if there is an img directory, if not then create one.
# =============================================================================

if os.path.isdir(dir+'/img') == True:
    pass
else:
    os.system('mkdir {}/img'.format(dir))

# =============================================================================
# Extract the snapshot filename 
# =============================================================================

for idx,lines in enumerate(os.listdir(dir)):
    if lines.find('snapshot') != -1:
        snapshot_file_name = os.listdir(dir)[idx]

# =============================================================================
# Load the data from the snapshot data file
# =============================================================================

with h5py.File('{}{}/snapshots.h5'.format(dir,snapshot_file_name), mode = 'r') as file:

    # =========================================================================
    # Velocity data
    # =========================================================================
    
    kinetic_spectrum = np.copy(file['tasks']["kinetic_spectrum"])[-snap_t:,:,:,:]
    u = np.copy(file['tasks']["u"])[-snap_t:,:,:,:]
    v = np.copy(file['tasks']["v"])[-snap_t:,:,:,:]
    w = np.copy(file['tasks']["w"])[-snap_t:,:,:,:]
    T = np.copy(file['tasks']["T"])[-snap_t:,:,:,:]

    # =========================================================================
    # Momentum X equationnp.sum(realFieldData, axis = 2)
    # =========================================================================
    
    x_pressure = np.copy(file['tasks']["x_pressure"])[-snap_t:,:,:,:]
    x_diffusion = np.copy(file['tasks']["x_diffusion"])[-snap_t:,:,:,:]
    x_coriolis = np.copy(file['tasks']["x_coriolis"])[-snap_t:,:,:,:]
    x_inertia = np.copy(file['tasks']["x_inertia"])[-snap_t:,:,:,:]

    # =========================================================================
    # Momentum Y equation
    # =========================================================================
    
    y_pressure = np.copy(file['tasks']["y_pressure"])[-snap_t:,:,:,:]
    y_diffusion = np.copy(file['tasks']["y_diffusion"])[-snap_t:,:,:,:]
    y_coriolis = np.copy(file['tasks']["y_coriolis"])[-snap_t:,:,:,:]
    y_inertia =	np.copy(file['tasks']["y_inertia"])[-snap_t:,:,:,:]

    # =========================================================================
    # Momentum Z equation
    # =========================================================================
    
    z_pressure = np.copy(file['tasks']["z_pressure"])[-snap_t:,:,:,:]
    z_diffusion = np.copy(file['tasks']["z_diffusion"])[-snap_t:,:,:,:]
    z_inertia = np.copy(file['tasks']["z_inertia"])[-snap_t:,:,:,:]
    z_buoyancy = np.copy(file['tasks']["z_bouyancy"])[-snap_t:,:,:,:]
    
    # =========================================================================
    # Vorticity X equation
    # =========================================================================
    
    vorticity_x_diffusion = np.copy(file['tasks']["vorticity_x_diffusion"])[-snap_t:,:,:,:]
    vorticity_x_coriolis = np.copy(file['tasks']["vorticity_x_coriolis"])[-snap_t:,:,:,:]
    vorticity_x_bouyancy = np.copy(file['tasks']["vorticity_x_bouyancy"])[-snap_t:,:,:,:]
    vorticity_x_inertia = np.copy(file['tasks']["vorticity_x_inertia"])[-snap_t:,:,:,:]

    # =========================================================================
    # Vorticity Y equation
    # =========================================================================
    
    vorticity_y_diffusion = np.copy(file['tasks']["vorticity_y_diffusion"])[-snap_t:,:,:,:]
    vorticity_y_coriolis = np.copy(file['tasks']["vorticity_y_coriolis"])[-snap_t:,:,:,:]
    vorticity_y_bouyancy = np.copy(file['tasks']["vorticity_y_bouyancy"])[-snap_t:,:,:,:]
    vorticity_y_inertia =	np.copy(file['tasks']["vorticity_y_inertia"])[-snap_t:,:,:,:]

    # =========================================================================
    # Vorticity Z equation
    # =========================================================================
    
    vorticity_z_diffusion = np.copy(file['tasks']["vorticity_z_diffusion"])[-snap_t:,:,:,:]
    vorticity_z_inertia = np.copy(file['tasks']["vorticity_z_inertia"])[-snap_t:,:,:,:]
    vorticity_z_coriolis = np.copy(file['tasks']["vorticity_z_coriolis"])[-snap_t:,:,:,:]
    

# =============================================================================
# All my functions which I use through out the script
# =============================================================================

def plotASlice(Field):
    """ 
    Parameters
    ----------
    Field : Dedalus Field object

    Raises
    ------
    Must be working with 3D data.

    Returns
    -------
    None.
    
    Notes
    -----
    Creates several plots looking at the length scales in the given field.
    """

    realFieldData = Field['g']
    fieldName = Field.name
    midPoint = np.shape(realFieldData)[-1] // 2

    if len(np.shape(realFieldData)) != 3:
        raise Exception('Expecting 3d data.')
    
    slice = np.sum(realFieldData, axis = 2)
    sliceXAverage = np.average(slice, axis = 0)
    sliceYAverage = np.average(slice, axis = 1)
    horizontalDomain = np.linspace(0, 2, len(sliceYAverage)) 

    fig,ax = plt.subplots(2,2, figsize=(10,10))
    ax[0][0].imshow(slice, cmap = 'coolwarm')
    ax[0][0].set_title('Average {} Field'.format(fieldName))
    ax[1][0].plot(horizontalDomain, sliceXAverage, lw = spectrumlw, color = CB91_Blue)
    ax[1][0].set_title('X Average')   
    ax[0][1].plot(horizontalDomain, sliceYAverage, lw = spectrumlw, color = CB91_Violet)
    ax[0][1].set_title('Y Average')
    spectraX = (fft(sliceXAverage).real**2 + fft(sliceXAverage).imag**2)**0.5
    spectraY = (fft(sliceYAverage).real**2 + fft(sliceYAverage).imag**2)**0.5
    midPoint = len(spectraX) // 2
    ax[1][1].plot(spectraX[:midPoint], lw = spectrumlw, color = CB91_Blue)
    ax[1][1].plot(spectraY[:midPoint], lw = spectrumlw, color = CB91_Violet) 
    ax[1][1].set_title('Spectra')
    plt.savefig('{}/img/plotASlice.eps'.format(dir), dpi=500)
    plt.show()
    
def determineRoot(derivative,z):
        """ 
         

        Parameters
        ----------
        derivative : The derivative of a profile.
        z : The z domain.

        Returns
        -------
        upper : Upper boundarie.
        lower : Lower boundarie.
        avg : Average of the two.
        avg_points : Average number of points in both.

        """

        upper,lower,avg = 0,0,0
        zero_crossings = np.where(np.diff(np.sign(derivative)))[0]

        # =========================================================================
        # Determine the lower crossing 
        # =========================================================================
        
        x,y = [z[zero_crossings[0]],z[zero_crossings[0]+1]], [derivative[zero_crossings[0]], derivative[zero_crossings[0]+1]]
        m,b = np.polyfit(y,x,1)
        lower = b
        
        # =========================================================================
        # Determine the upper crossing 
        # =========================================================================
        
        x,y = [z[zero_crossings[-1]],z[zero_crossings[-1]+1]], [derivative[zero_crossings[-1]],derivative[zero_crossings[-1]+1]]
        m,b = np.polyfit(y,x,1)
        upper = b
        avg = (lower + (1-upper))/2

        # =========================================================================
        # Calculate the average number of points in the boundaries.
        # =========================================================================
        
        avg_points = zero_crossings[0]*mask

        return upper, lower, avg, avg_points

def derivative(data,z):
    """ 
    

    Parameters
    ----------
    data : 1d data.
    z : z domain array.

    Returns
    -------
    Derivative with respect to z.

    """
    return np.gradient(data,z)

def plotMidPlane(Field):
    """ 
    Parameters
    ----------
    Field : Dedalus Field object

    Raises
    ------
    Must be working with 3D data.

    Returns
    -------
    None.
    
    Notes
    -----
    Created a imshow plot of the mid section of the domain.
    """

    realFieldData = Field['g']
    fieldName = Field.name
    midPoint = np.shape(realFieldData)[-1] // 2

    if len(np.shape(realFieldData)) != 3:
        raise Exception('Expecting 3d data.')
    
    rotatedData = np.rot90(realFieldData, k=1, axes = (2,0))
    midPointSlice = rotatedData[midPoint]
    fig = plt.figure(figsize = (6,6))
    plt.xticks([])
    plt.yticks([])
    plt.imshow(midPointSlice, cmap = 'bwr', interpolation = 'gaussian')
    plt.show()
    
def computeRMS(Fx, Fy, Fz):
    """ 
    

    Parameters
    ----------
    Fx : FORCE IN X.
    Fy : FORCE IN Y.
    Fz : FORCE IN Z.

    Returns
    -------
    Removes any mean profile from the flow and returns the RMS of this 
    new, meanless, field.

    """

    return np.sqrt(Fx**2 + Fy**2 + Fz**2)


def removeMeanProfile(Force):
    """ 
    

    Parameters
    ----------
    Force : FORCE.

    Returns
    -------
    Removes any mean profile found in the flow (Force).
    """
    
    rotatedForce = np.rot90(Force, k = 1, axes = (2,0))
    meanProfileOfForce = [np.average(planes) for planes in rotatedForce]
    meanProfile3DArray = np.ones_like(rotatedForce)
    
    for idx,plane in enumerate(meanProfile3DArray):
        meanProfile3DArray[idx] = meanProfile3DArray[idx]*meanProfileOfForce[idx]

    return np.rot90(rotatedForce - meanProfile3DArray, k = -1, axes = (2,0)) 

def horizontalAverage(ForceRMS):
    """ 
    

    Parameters
    ----------
    ForceRMS : Root mean square of the force.
ffns
    -------yVorticityViscositySpectrum,
    The averaged of each force over x and y, leaving a 1d profile against z. 
    As done in (Guzman, 2021).

    """
    profile = []
    rotatedForce = np.rot90(ForceRMS, k=1, axes=(2,0))
 
    for lines in rotatedForce:
        profile.append(np.average(lines))

    return np.array(profile)


def removeBoundaries(Mask, Force):
    """ 
    

    Parameters
    ----------
    Mask : An array of ones and zeroes which removes the boundaries.
    Force : The force from which we want to remove the boundaries.

    Returns
    -------
    None: preformed in place.

    """

    RotatedForce = np.rot90(Force['g'], k=1, axes = (2,0))
    MaskedForce = RotatedForce*Mask
    Force['g'] = np.rot90(MaskedForce,	k=-1, axes = (2,0))

# =============================================================================
# Create a dedalus domain which is the same as the domain of the problem so
# we can use in-built dedalus functions to preform integration and compute
# spectra.
# =============================================================================

Nx = N 
Ny = N
Lx = Ly = 1
Lz = 1
Nz = int(N/Lx)

x_basis = de.Fourier('x', Nx, interval = (0,Lx), dealias=3/2)
y_basis = de.Fourier('y', Ny, interval = (0,Ly), dealias=3/2)
z_basis = de.Chebyshev('z', Nz, interval = (-Lz/2,Lz/2), dealias =3/2)
domain = de.Domain([x_basis, y_basis, z_basis], grid_dtype=np.float64)

# =============================================================================
# Create fields in dedalus so we can use dedalus tools for computations, 
# retain spectral accuracy e.c.t.
# =============================================================================

# =============================================================================
# Momentum Equation
# =============================================================================

Viscosity = domain.new_field(name='viscosity')
Coriolis = domain.new_field(name='coriolis')
ACoriolis = domain.new_field(name='Acoriolis')
Inertia = domain.new_field(name='inertia')
Buoyancy = domain.new_field(name='buoyancy')
Pressure = domain.new_field(name='pressure')


# =============================================================================
# Avg the horizontal velocity profile over time
# Create an array which contains only the last N points, then avg over this.
# =============================================================================

avg_u_h = np.average(np.array(U_h[-transient:,0,0,:]), axis=0) 
upper_viscous_boundary, lower_viscous_boundary, avg_viscous_boundary, avg_points = determineRoot(derivative(avg_u_h,z),z)

# =============================================================================
# Construct a matrix which is 1 everywhere other than in the viscous boundary 
# where it is 0. This can be used to remove the viscous boundaries from my 
# force plots.
# =============================================================================

Mask = np.rot90(np.ones(np.shape(x_diffusion[-1]),dtype=np.int32), k=1 , axes=(2,0))

for idx,slices in enumerate(Mask):
    if idx < avg_points or idx >= len(z) - avg_points:
        Mask[idx] = np.zeros(np.shape(slices))

# =============================================================================
# Arrays containing the time series for the horizontal average
# =============================================================================

horizontalAvgViscosityTimeSeries = []
horizontalAvgCoriolisTimeSeries = []
horizontalAvgBuoyancyTimeSeries = []
horizontalAvgInertiaTimeSeries = []
horizontalAvgPressureTimeSeries = []
horizontalAvgACoriolisTimeSeries = []
blankMatrix = np.zeros(np.shape(z_diffusion[-1]), dtype=np.int32)

for idx in range(1,snap_t+1):

    # =============================================================================
    # Load data
    # =============================================================================
    
    Viscosity['g'] = computeRMS(x_diffusion[-idx], y_diffusion[-idx], z_diffusion[-idx])
    Coriolis['g'] = computeRMS(x_coriolis[-idx], y_coriolis[-idx], blankMatrix)
    Inertia['g'] = computeRMS(x_inertia[-idx], y_inertia[-idx], z_inertia[-idx])
    Buoyancy['g'] = computeRMS(blankMatrix, blankMatrix, removeMeanProfile(z_buoyancy[-idx]))
    Pressure['g'] = computeRMS(removeMeanProfile(x_pressure[-idx]), removeMeanProfile(y_pressure[-idx]), removeMeanProfile(z_pressure[-idx]))
    ACoriolis['g'] = computeRMS(x_pressure[-idx] + x_coriolis[-idx], y_pressure[-idx] + y_coriolis[-idx], z_pressure[-idx] - blankMatrix)

    # =========================================================================
    # Remove the boundaries
    # =========================================================================
    
    removeBoundaries(Mask, Viscosity)
    removeBoundaries(Mask, Pressure)
    removeBoundaries(Mask, Coriolis)
    removeBoundaries(Mask, Inertia)
    removeBoundaries(Mask, ACoriolis)
    removeBoundaries(Mask, Buoyancy)

    # =========================================================================
    # Compute the horizontal average
    # =========================================================================
        
    horizontalAvgViscosityTimeSeries.append(abs(horizontalAverage(Viscosity['g'])))
    horizontalAvgCoriolisTimeSeries.append(abs(horizontalAverage(Coriolis['g'])))
    horizontalAvgInertiaTimeSeries.append(abs(horizontalAverage(Inertia['g'])))
    horizontalAvgBuoyancyTimeSeries.append(abs(horizontalAverage(Buoyancy['g'])))
    horizontalAvgPressureTimeSeries.append(abs(horizontalAverage(Pressure['g'])))
    horizontalAvgACoriolisTimeSeries.append(abs(horizontalAverage(ACoriolis['g'])))

    print('Coriolis: {:.2e}, Pressure: {:.2e}, Viscosity: {:.2e}, Inertia: {:.2e}, Buoyancy: {:.2e}, Ageostrophic: {:.2e}'.format(np.sum(Coriolis['g']), np.sum(Pressure['g']),np.sum(Viscosity['g']), np.sum(Inertia['g']), np.sum(Buoyancy['g']), np.sum(ACoriolis['g'])))

# =============================================================================
# Time avg the profiles
# =============================================================================

ViscosityProfile = np.average(np.array(horizontalAvgViscosityTimeSeries), axis=0)
InertiaProfile = np.average(np.array(horizontalAvgInertiaTimeSeries), axis=0)
BuoyancyProfile = np.average(np.array(horizontalAvgBuoyancyTimeSeries), axis=0)
CoriolisProfile = np.average(np.array(horizontalAvgCoriolisTimeSeries), axis=0)
PressureProfile = np.average(np.array(horizontalAvgPressureTimeSeries), axis=0)
ACoriolisProfile = np.average(np.array(horizontalAvgACoriolisTimeSeries), axis=0)

# =============================================================================
# Plot the horizontal profiles 
# =============================================================================

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

# =============================================================================
# ## Plot mid height values like (Guzman, 2021) ##
# ============================================================================= 

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

print('-'*120)
print('Viscosity: {:.3e}, Coriolis: {:.3e}, Inertia: {:.3e}, Buoyancy: {:.3e}, Pressure: {:.3e}.'.format(ViscosityProfile[mid_point], CoriolisProfile[mid_point], InertiaProfile[mid_point], BuoyancyProfile[mid_point], PressureProfile[mid_point]))
print('-'*120)
