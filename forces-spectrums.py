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

def findIntersection(arrOne, arrTwo):
    """
    

    Parameters
    ----------
    arrOne : Spectra one.
    arrTwo : Spectra two.

    Returns
    -------
    Returns the intersection between two lines by considering the difference 
    between the two and returning where this line crosses 0. If statement 
    catches the case where the lines dont intersect.

    """
    
    difference = arrOne - arrTwo
    minimumDifferenceIndex = np.argmin(difference)
    
    signChange = np.where(np.diff(np.sign(difference)) != 0)[0] + 1
    
    if len(signChange) > 0:
        m,b = np.polyfit([signChange[0] - 1, signChange[0]] , [difference[signChange[0] - 1], difference[signChange[0]]], 1)
        return -b/m, abs(b)

    return 1,1

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
    
def compareSpectrum(Field):
    """
    

    Parameters
    ----------
    Field : TA dedalus field object

    Returns
    -------
    None.
    
    Notes
    -----
    Integrates the field in real space and then computes a summation of the 
    spectrum to show that they are indeed the same.

    """
    
    integratedRealField = Field.integrate()['g'][0][0][0]
    spectrum = Field['c']
    runningSum = 0
    realZerothSpectrum = spectrum[0][0][:].real
 
    for idx,lines in enumerate(realZerothSpectrum):
         
        if idx == 0:
            runningSum += 2*spectrum[0][0][idx]

        elif idx % 2 == 0:
            runningSum -= 2*spectrum[0][0][idx] / (idx**2 - 1)  
        else:
           pass

    print(2*runningSum, integratedRealField)

def computeRMS(Fx, Fy, Fz):
    """ 
    

    Parameters
    ----------
    Fx : FORCE IN X.
    Fy : FORCE IN Y.
    Fz : FORCE IN Z.

    ReturnsyVorticityViscositySpectrum,
    -------
    RMS of the force. Computes: sqrt(Fx^2 + Fy^2 + Fz^2) as used in (Andres,2021).

    """

    return Fx + Fy + Fz#np.sqrt(Fx**2 + Fy**2 + Fz**2)


def computeSpectrum(Force):
    """ 
    

    Parameters
    ----------
    Force : Field object from dedalus.

    Returns
    -------
    The spectrum of the force averaged over two directions.

    """
    
    ForceSpectrum = (Force['c'].imag**2 +  Force['c'].real**2)**0.5
    z_avg = np.sum(ForceSpectrum, axis = 2) # average over z 
    return np.sum(z_avg, axis = 1) 

def horizontalAverage(ForceRMS):
    """ 
    

    Parameters
    ----------
    ForceRMS : Root mean square of the force.

    Returns
    -------yVorticityViscositySpectrum,
    The averaged of each force over x and y, leaving a 1d profile against z. 
    As done in (Guzman, 2021).

    """
    profile = []
    rotatedForce = np.rot90(ForceRMS, k=1, axes=(2,0))
 
    for lines in rotatedForce:
        profile.append(np.sum(lines))

    return np.array(profile)

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
Lx = Ly = 2
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
# Full Vorticity Equation
# =============================================================================

vorticityViscosity = domain.new_field(name='vorticityViscosity')
vorticityCoriolis = domain.new_field(name='vorticityCoriolis')
vorticityInertia = domain.new_field(name='vorticityInertia')
vorticityBuoyancy = domain.new_field(name='vorticityBuoyancy')

# =============================================================================
# X Vorticity Equation
# =============================================================================

xVorticityViscosity = domain.new_field(name='xVorticityViscosity')
xVorticityCoriolis = domain.new_field(name='xVorticityCoriolis')
xVorticityInertia = domain.new_field(name='xVorticityInertia')
xVorticityBuoyancy = domain.new_field(name='xVorticityBuoyancy')

# =============================================================================
# Y Vorticity Equation
# =============================================================================

yVorticityViscosity = domain.new_field(name='yVorticityViscosity')
yVorticityCoriolis = domain.new_field(name='yVorticityCoriolis')
yVorticityInertia = domain.new_field(name='yVorticityInertia')
yVorticityBuoyancy = domain.new_field(name='yVorticityBuoyancy')

# =============================================================================
# X and Y averaged Vorticity Equation
# =============================================================================

xyVorticityViscosity = domain.new_field(name='xyVorticityViscosity')
xyVorticityCoriolis = domain.new_field(name='xyVorticityCoriolis')
xyVorticityInertia = domain.new_field(name='xyVorticityInertia')
xyVorticityBuoyancy = domain.new_field(name='xyVorticityBuoyancy')

# =============================================================================
# Velocity fields
# =============================================================================

Kinetic = domain.new_field(name='kinetic')
U = domain.new_field(name='u')
V = domain.new_field(name='v')
W = domain.new_field(name='w')
Temperature = domain.new_field(name='Temperature')
Temperature['g'] = T[-1]

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
# Arrays containing the time series for the spectra
# =============================================================================

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

vorticityViscosityTimeSeries = []
vorticityCoriolisTimeSeries = []
vorticityBuoyancyTimeSeries = []
vorticityInertiaTimeSeries = []

xVorticityViscosityTimeSeries = []
xVorticityCoriolisTimeSeries = []
xVorticityBuoyancyTimeSeries = []
xVorticityInertiaTimeSeries = []

yVorticityViscosityTimeSeries = []
yVorticityCoriolisTimeSeries = []
yVorticityBuoyancyTimeSeries = []
yVorticityInertiaTimeSeries = []

xyVorticityViscosityTimeSeries = []
xyVorticityCoriolisTimeSeries = []
xyVorticityBuoyancyTimeSeries = []
xyVorticityInertiaTimeSeries = []

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
    Buoyancy['g'] = computeRMS(blankMatrix, blankMatrix, z_buoyancy[-idx])
    Pressure['g'] = computeRMS(x_pressure[-idx], y_pressure[-idx], z_pressure[-idx])
    Kinetic['c'] = kinetic_spectrum[-idx]
    ACoriolis['g'] = computeRMS(x_pressure[-idx] + x_coriolis[-idx], y_pressure[-idx] + y_coriolis[-idx], z_pressure[-idx] - blankMatrix)
    U['g'] = u[-idx]
    V['g'] = v[-idx]
    W['g'] = w[-idx]
    
    vorticityViscosity['g'] = computeRMS(vorticity_x_diffusion[-idx], vorticity_y_diffusion[-idx], vorticity_z_diffusion[-idx])
    vorticityCoriolis['g'] = computeRMS(vorticity_x_coriolis[-idx], vorticity_y_coriolis[-idx], vorticity_z_coriolis[-idx])
    vorticityInertia['g'] = computeRMS(vorticity_x_inertia[-idx], vorticity_y_inertia[-idx], vorticity_z_inertia[-idx])
    vorticityBuoyancy['g'] = computeRMS(vorticity_x_bouyancy[-idx], vorticity_y_bouyancy[-idx], blankMatrix)
    
    xVorticityViscosity['g'] = vorticity_x_diffusion[-idx]
    xVorticityCoriolis['g'] = vorticity_x_coriolis[-idx]
    xVorticityInertia['g'] = vorticity_x_inertia[-idx]
    xVorticityBuoyancy['g'] = vorticity_x_bouyancy[-idx]  
 
    yVorticityViscosity['g'] = vorticity_y_diffusion[-idx]
    yVorticityCoriolis['g'] = vorticity_y_coriolis[-idx]
    yVorticityInertia['g'] = vorticity_y_inertia[-idx]
    yVorticityBuoyancy['g'] = vorticity_y_bouyancy[-idx]
 
    xyVorticityViscosity['g'] = computeRMS(vorticity_x_diffusion[-idx], vorticity_y_diffusion[-idx], blankMatrix)
    xyVorticityCoriolis['g'] = computeRMS(vorticity_x_coriolis[-idx], vorticity_y_coriolis[-idx], blankMatrix)
    xyVorticityInertia['g'] = computeRMS(vorticity_x_inertia[-idx], vorticity_y_inertia[-idx], blankMatrix)
    xyVorticityBuoyancy['g'] = computeRMS(vorticity_x_bouyancy[-idx], vorticity_y_bouyancy[-idx], blankMatrix)

    # =========================================================================
    # Remove the boundaries
    # =========================================================================
    
    removeBoundaries(Mask, Viscosity)
    removeBoundaries(Mask, Pressure)
    removeBoundaries(Mask, Coriolis)
    removeBoundaries(Mask, Inertia)
    removeBoundaries(Mask, ACoriolis)
    removeBoundaries(Mask, Buoyancy)
    removeBoundaries(Mask, Kinetic)
    removeBoundaries(Mask, U)
    removeBoundaries(Mask, V)
    removeBoundaries(Mask, W)
    removeBoundaries(Mask, vorticityViscosity)
    removeBoundaries(Mask, vorticityCoriolis)
    removeBoundaries(Mask, vorticityInertia)
    removeBoundaries(Mask, vorticityBuoyancy)
    removeBoundaries(Mask, xVorticityViscosity)
    removeBoundaries(Mask, xVorticityCoriolis)
    removeBoundaries(Mask, xVorticityInertia)
    removeBoundaries(Mask, xVorticityBuoyancy)
    removeBoundaries(Mask, yVorticityViscosity)
    removeBoundaries(Mask, yVorticityCoriolis)
    removeBoundaries(Mask, yVorticityInertia)
    removeBoundaries(Mask, yVorticityBuoyancy)
    removeBoundaries(Mask, xyVorticityViscosity)
    removeBoundaries(Mask, xyVorticityCoriolis)
    removeBoundaries(Mask, xyVorticityInertia)
    removeBoundaries(Mask, xyVorticityBuoyancy)

    
    # =========================================================================
    # Compute the spectrum
    # =========================================================================
    
    ViscosityTimeSeries.append(computeSpectrum(Viscosity))
    CoriolisTimeSeries.append(computeSpectrum(Coriolis))
    InertiaTimeSeries.append(computeSpectrum(Inertia))
    BuoyancyTimeSeries.append(computeSpectrum(Buoyancy))
    PressureTimeSeries.append(computeSpectrum(Pressure))
    KineticTimeSeries.append(computeSpectrum(Kinetic))
    UTimeSeries.append(computeSpectrum(U))
    VTimeSeries.append(computeSpectrum(V))
    WTimeSeries.append(computeSpectrum(W))
    ACoriolisTimeSeries.append(computeSpectrum(ACoriolis))
    vorticityViscosityTimeSeries.append(computeSpectrum(vorticityViscosity))
    vorticityCoriolisTimeSeries.append(computeSpectrum(vorticityCoriolis))
    vorticityInertiaTimeSeries.append(computeSpectrum(vorticityInertia))
    vorticityBuoyancyTimeSeries.append(computeSpectrum(vorticityBuoyancy))

    xVorticityViscosityTimeSeries.append(computeSpectrum(xVorticityViscosity))
    xVorticityCoriolisTimeSeries.append(computeSpectrum(xVorticityCoriolis))
    xVorticityInertiaTimeSeries.append(computeSpectrum(xVorticityInertia))
    xVorticityBuoyancyTimeSeries.append(computeSpectrum(xVorticityBuoyancy))

    yVorticityViscosityTimeSeries.append(computeSpectrum(yVorticityViscosity))
    yVorticityCoriolisTimeSeries.append(computeSpectrum(yVorticityCoriolis))
    yVorticityInertiaTimeSeries.append(computeSpectrum(yVorticityInertia))
    yVorticityBuoyancyTimeSeries.append(computeSpectrum(yVorticityBuoyancy))

    xyVorticityViscosityTimeSeries.append(computeSpectrum(xyVorticityViscosity))
    xyVorticityCoriolisTimeSeries.append(computeSpectrum(xyVorticityCoriolis))
    xyVorticityInertiaTimeSeries.append(computeSpectrum(xyVorticityInertia))
    xyVorticityBuoyancyTimeSeries.append(computeSpectrum(xyVorticityBuoyancy))

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

compareSpectrum(Viscosity)

# =============================================================================
# Time avg the spectrums 
# =============================================================================

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

vorticityViscositySpectrum = np.average(np.array(vorticityViscosityTimeSeries), axis=0)
vorticityInertiaSpectrum = np.average(np.array(vorticityInertiaTimeSeries), axis=0)
vorticityBuoyancySpectrum = np.average(np.array(vorticityBuoyancyTimeSeries), axis=0)
vorticityCoriolisSpectrum = np.average(np.array(vorticityCoriolisTimeSeries), axis=0)

xVorticityViscositySpectrum = np.average(np.array(xVorticityViscosityTimeSeries), axis=0)
xVorticityInertiaSpectrum = np.average(np.array(xVorticityInertiaTimeSeries), axis=0)
xVorticityBuoyancySpectrum = np.average(np.array(xVorticityBuoyancyTimeSeries), axis=0)
xVorticityCoriolisSpectrum = np.average(np.array(xVorticityCoriolisTimeSeries), axis=0)

yVorticityViscositySpectrum = np.average(np.array(yVorticityViscosityTimeSeries), axis=0)
yVorticityInertiaSpectrum = np.average(np.array(yVorticityInertiaTimeSeries), axis=0)
yVorticityBuoyancySpectrum = np.average(np.array(yVorticityBuoyancyTimeSeries), axis=0)
yVorticityCoriolisSpectrum = np.average(np.array(yVorticityCoriolisTimeSeries), axis=0)

xyVorticityViscositySpectrum = np.average(np.array(xyVorticityViscosityTimeSeries), axis=0)
xyVorticityInertiaSpectrum = np.average(np.array(xyVorticityInertiaTimeSeries), axis=0)
xyVorticityBuoyancySpectrum = np.average(np.array(xyVorticityBuoyancyTimeSeries), axis=0)
xyVorticityCoriolisSpectrum = np.average(np.array(xyVorticityCoriolisTimeSeries), axis=0)

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
# Plotting the results
# =============================================================================

plotASlice(W)

plotMidPlane(Temperature)

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
plt.ylim( min(yVorticityCoriolisSpectrum) // 10, max(yVorticityCoriolisSpectrum)*3)
viscosityIndex, viscosityIntersection = findIntersection(ViscositySpectrum, ACoriolisSpectrum)
print('The intersection between viscosity and the Ageostrophic Coriolis force in the force spectra is: {:.2f}.'.format(viscosityIndex))
plt.plot(np.ones(100)*viscosityIndex, np.linspace(0, viscosityIntersection, 100), 'k--')
plt.annotate('$L_v = ${:.1f}'.format(viscosityIndex), (viscosityIndex // 2, max(ViscositySpectrum) // 10))
inertiaIndex, inertiaIntersection = findIntersection(InertiaSpectrum, ACoriolisSpectrum)
plt.plot(np.ones(100)*inertiaIndex, np.linspace(0, inertiaIntersection, 100), 'k--')
plt.annotate('$L_I = ${:.1f}'.format(inertiaIndex), (inertiaIndex // 2, max(ViscositySpectrum) // 20))
legend_properties = {'weight':'bold'}
plt.title('Force Spectra')
plt.legend(ncol=2, fontsize=14,prop=legend_properties,frameon=False)
plt.savefig('{}/img/ForceSpectrum.eps'.format(dir), dpi=500)
plt.show()

fig = plt.figure(figsize=(10,10))
plt.plot(range(len(ViscositySpectrum)), vorticityViscositySpectrum, label = '$\\omega_v$', color = ViscosityColour, lw=spectrumlw)
plt.plot(range(len(ViscositySpectrum)), vorticityCoriolisSpectrum, label = '$\\omega_C$', color = CoriolisColour, lw=spectrumlw)
plt.plot(range(len(InertiaSpectrum)), vorticityInertiaSpectrum, label = '$\\omega_I$', color = InertiaColour, lw=spectrumlw)
plt.plot(range(len(BuoyancySpectrum)), vorticityBuoyancySpectrum, label = '$\\omega_B$', color = BuoyancyColour, lw=spectrumlw)
plt.xscale("log")
plt.yscale("log")
plt.xlabel('$K_z$')
plt.title('Vorticity Equation')
plt.ylabel('Magnitude')
plt.xlim(1,64)
plt.ylim( min(vorticityCoriolisSpectrum) // 10, max(vorticityCoriolisSpectrum)*3)
viscosityIndex, viscosityIntersection = findIntersection(vorticityViscositySpectrum, vorticityCoriolisSpectrum)
print('The intersection between viscosity and the Ageostrophic Coriolis force in the force spectra is: {:.2f}.'.format(viscosityIndex))
plt.plot(np.ones(100)*viscosityIndex, np.linspace(0, viscosityIntersection, 100), 'k--')
plt.annotate('$L_v = ${:.1f}'.format(viscosityIndex), (viscosityIndex // 2, max(ViscositySpectrum) // 10))
inertiaIndex, inertiaIntersection = findIntersection(vorticityInertiaSpectrum, vorticityCoriolisSpectrum)
plt.plot(np.ones(100)*inertiaIndex, np.linspace(0, inertiaIntersection, 100), 'k--')
plt.annotate('$L_I = ${:.1f}'.format(inertiaIndex), (inertiaIndex // 2, max(ViscositySpectrum) // 20))
plt.legend(ncol=2, fontsize=14,prop=legend_properties,frameon=False)
plt.savefig('{}/img/vorticityForceSpectrum.eps'.format(dir), dpi=500)
plt.show()

fig = plt.figure(figsize=(10,10))
plt.plot(range(len(ViscositySpectrum)), yVorticityViscositySpectrum, label = '$\\omega_v$', color = ViscosityColour, lw=spectrumlw)
plt.plot(range(len(ViscositySpectrum)), yVorticityCoriolisSpectrum, label = '$\\omega_C$', color = CoriolisColour, lw=spectrumlw)
plt.plot(range(len(InertiaSpectrum)), yVorticityInertiaSpectrum, label = '$\\omega_I$', color = InertiaColour, lw=spectrumlw)
plt.plot(range(len(BuoyancySpectrum)), yVorticityBuoyancySpectrum, label = '$\\omega_B$', color = BuoyancyColour, lw=spectrumlw)
plt.xscale("log")
plt.title('Y - Vorticity Equation')
plt.yscale("log")
plt.xlabel('$K_z$')
plt.ylabel('Magnitude')
viscosityIndex, viscosityIntersection = findIntersection(yVorticityViscositySpectrum, yVorticityCoriolisSpectrum)
print('The intersection between viscosity and the Ageostrophic Coriolis force in the force spectra is: {:.2f}.'.format(viscosityIndex))
plt.plot(np.ones(100)*viscosityIndex, np.linspace(0, viscosityIntersection, 100), 'k--')
plt.annotate('$L_v = ${:.1f}'.format(viscosityIndex), (viscosityIndex // 2, max(ViscositySpectrum) // 10))
inertiaIndex, inertiaIntersection = findIntersection(yVorticityInertiaSpectrum, yVorticityCoriolisSpectrum)
plt.plot(np.ones(100)*inertiaIndex, np.linspace(0, inertiaIntersection, 100), 'k--')
plt.annotate('$L_I = ${:.1f}'.format(inertiaIndex), (inertiaIndex // 2, max(ViscositySpectrum) // 20))
plt.xlim(1,64) 
plt.ylim( min(yVorticityCoriolisSpectrum) // 10, max(yVorticityCoriolisSpectrum)*3)
plt.legend(ncol=2, fontsize=14,prop=legend_properties,frameon=False)
plt.savefig('{}/img/yVorticityForceSpectrum.eps'.format(dir), dpi=500)
plt.show()

fig = plt.figure(figsize=(10,10))
plt.plot(range(len(ViscositySpectrum)), xyVorticityViscositySpectrum, label = '$\\omega_v$', color = ViscosityColour, lw=spectrumlw)
plt.plot(range(len(ViscositySpectrum)), xyVorticityCoriolisSpectrum, label = '$\\omega_C$', color = CoriolisColour, lw=spectrumlw)
plt.plot(range(len(InertiaSpectrum)), xyVorticityInertiaSpectrum, label = '$\\omega_I$', color = InertiaColour, lw=spectrumlw)
plt.plot(range(len(BuoyancySpectrum)), xyVorticityBuoyancySpectrum, label = '$\\omega_B$', color = BuoyancyColour, lw=spectrumlw)
plt.xscale("log")
plt.title('X - Y Averaged Vorticity Equation')
plt.yscale("log")
plt.xlabel('$K_z$')
plt.ylabel('Magnitude')
plt.xlim(1,64)
index,intersection = findIntersection(xyVorticityViscositySpectrum, xyVorticityCoriolisSpectrum)
print('The intersection between viscosity and the Coriolis force in the X-Y avg vorticity equation is: {:.2f}.'.format(index))
plt.annotate('$L_\perp = ${:.1f}'.format(index), (index // 2, max(xyVorticityViscositySpectrum) // 10))
plt.plot(np.ones(100)*index, np.linspace(0, intersection, 100), 'k--')
plt.legend(ncol=2, fontsize=14,prop=legend_properties,frameon=False)
plt.savefig('{}/img/xyVorticityForceSpectrum.eps'.format(dir), dpi=500)
plt.show()

fig = plt.figure(figsize=(10,10))
plt.plot(range(len(ViscositySpectrum)), xVorticityViscositySpectrum, label = '$\\omega_v$', color = ViscosityColour, lw=spectrumlw)
plt.plot(range(len(ViscositySpectrum)), xVorticityCoriolisSpectrum, label = '$\\omega_C$', color = CoriolisColour, lw=spectrumlw)
plt.plot(range(len(InertiaSpectrum)), xVorticityInertiaSpectrum, label = '$\\omega_I$', color = InertiaColour, lw=spectrumlw)
plt.plot(range(len(BuoyancySpectrum)), xVorticityBuoyancySpectrum, label = '$\\omega_B$', color = BuoyancyColour, lw=spectrumlw)
plt.xscale("log")
plt.title('X - Vorticity Equation')
plt.yscale("log")
plt.xlabel('$K_z$')
plt.ylabel('Magnitude')
plt.xlim(1,64)
index,intersection = findIntersection(xVorticityViscositySpectrum, xVorticityCoriolisSpectrum)
print('The intersection between viscosity and the Coriolis force in the X vorticity equation is: {:.2f}.'.format(index))
plt.plot(np.ones(100)*index, np.linspace(0, intersection, 100), 'k--')
plt.legend(ncol=2, fontsize=14,prop=legend_properties,frameon=False)
plt.annotate('$L_\perp = ${:.1f}'.format(index), (index // 2, max(xVorticityViscositySpectrum) // 10))
plt.savefig('{}/img/xVorticityForceSpectrum.eps'.format(dir), dpi=500)
plt.show()

print('The max wavenumber for k.e: {}, u: {}, v: {}, w: {}.'.format(np.argmax(KineticSpectrum[:64]), np.argmax(USpectrum[:64]), np.argmax(VSpectrum[:64]), np.argmax(WSpectrum[:64])))

# =============================================================================
# fig = plt.figure(figsize=(12,6))
# plt.plot(range(len(ViscositySpectrum)), KineticSpectrum, lw = 2, color = ViscosityColour, label = 'K.e')
# plt.plot(range(len(ViscositySpectrum)), USpectrum, lw = 2, color = CoriolisColour, label = 'u')
# plt.plot(range(len(ViscositySpectrum)), VSpectrum, lw = 2, color = PressureColour, label = 'v')
# plt.plot(range(len(ViscositySpectrum)), WSpectrum, lw = 2, color = ACColour, label = 'w')
# plt.xscale("log")
# plt.yscale("log")
# plt.xlabel('$K_z$')
# plt.ylabel('Magnitude')
# plt.xlim(1,64)
# legend_properties = {'weight':'bold'}
# plt.legend(ncol=2, fontsize=14,prop=legend_properties,frameon=False)
# plt.savefig('{}/img/KineticSpectrum.eps'.format(dir), dpi=500)
# plt.show()
# =============================================================================

# =============================================================================
# Plot the horizontal profiles 
# =============================================================================


# =============================================================================
# fig = plt.figure(figsize=(10,10))
# plt.plot(ViscosityProfile, z, color = ViscosityColour, label = '$F_V$', lw = 2)
# plt.plot(CoriolisProfile, z, color = CoriolisColour, label = '$F_C$', lw = 2)
# plt.plot(InertiaProfile, z, color = InertiaColour, label = '$F_I$', lw = 2)
# plt.plot(ACoriolisProfile, z, color = ACColour, label = '$F_{AC}$', lw = 2)
# plt.plot(BuoyancyProfile, z, color = BuoyancyColour, label = '$F_B$', lw = 2)
# plt.plot(PressureProfile, z, color = PressureColour, label = '$F_P$', lw = 2)
# legend_properties = {'weight':'bold'}
# plt.legend(ncol=2, fontsize=14,prop=legend_properties,frameon=False)
# plt.xlabel('Magnitude')
# plt.ylabel('z')
# plt.xscale('log')
# plt.savefig('{}/img/ForceProfiles.eps'.format(dir), dpi=500)
# plt.show()
# =============================================================================

# =============================================================================
# ## Plot mid height values like (Guzman, 2021) ##
# 
# mid_point = int(len(z)/2)
# 
# fig = plt.figure(figsize=(12,6))
# plt.xscale("log")
# plt.yscale("log")
# plt.scatter([Ra], ViscosityProfile[mid_point], label = 'Viscosity')
# plt.scatter([Ra], CoriolisProfile[mid_point],  label = 'Coriolis')
# plt.scatter([Ra], InertiaProfile[mid_point], label = 'Inertia')
# plt.scatter([Ra], BuoyancyProfile[mid_point], label = 'Buoyancy')
# plt.scatter([Ra], PressureProfile[mid_point], label = 'Pressure')
# plt.legend()
# plt.show()
# =============================================================================
