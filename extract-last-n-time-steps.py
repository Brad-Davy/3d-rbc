"""
   A script which extracts the last time step from a h5 file and saves
   it to a new file named: 'last-time-step.h5'. Note that this script will 
   throw an error if a file named last-time-step.h5 already exists in the 
   target directory.

Usage:
    forces_spectrums.py [--dir=<directory>] [--n=<Number of time steps>]
    forces_spectrums.py -h | --help

Options:
    -h --help                           Display this help message
    --dir=<directory>                   Directory [default: results/Ra_1-40e+08_Ek_1-00e-05_Pr_7-0_N_128_q_1-2_k_48_enhanced/]
    --n=<Number of time steps>          Number of time steps to be saved [default: 5]
"""

from docopt import docopt
import numpy as np
import matplotlib.pyplot as plt
import h5py
import os

%matplotlib 
# =============================================================================
# Extract the docopt arguments 
# =============================================================================

args = docopt(__doc__)
dataDirectory = str(args['--dir'])
numberOfTimeSteps = int(args['--n'])

# =============================================================================
# Extract the snapshot filename 
# =============================================================================

for idx,lines in enumerate(os.listdir(dataDirectory)):
    if lines.find('snapshot') != -1:
        snapshot_file_name = os.listdir(dataDirectory)[idx]
        
        
myData = h5py.File('{}{}/snapshots.h5'.format(dataDirectory,snapshot_file_name, str(numberOfTimeSteps)), mode = 'r')

u = np.copy(myData['tasks']['u'])[-numberOfTimeSteps:][:][:][:]
v = np.copy(myData['tasks']['v'])[-numberOfTimeSteps:][:][:][:]
w = np.copy(myData['tasks']['w'])[-numberOfTimeSteps:][:][:][:]
p = np.copy(myData['tasks']['p'])[-numberOfTimeSteps:][:][:][:]
T = np.copy(myData['tasks']['T'])[-numberOfTimeSteps:][:][:][:]
Tz = np.copy(myData['tasks']['Tz'])[-numberOfTimeSteps:][:][:][:]
uz = np.copy(myData['tasks']['uz'])[-numberOfTimeSteps:][:][:][:]
vz = np.copy(myData['tasks']['vz'])[-numberOfTimeSteps:][:][:][:]
wz = np.copy(myData['tasks']['wz'])[-numberOfTimeSteps:][:][:][:]


# =========================================================================
# Velocity data
# =========================================================================

kinetic_spectrum = np.copy(myData['tasks']["kinetic_spectrum"])[-numberOfTimeSteps:,:,:,:]

# =========================================================================
# Momentum X equationnp.sum(realFieldData, axis = 2)
# =========================================================================

x_pressure = np.copy(myData['tasks']["x_pressure"])[-numberOfTimeSteps:,:,:,:]
x_diffusion = np.copy(myData['tasks']["x_diffusion"])[-numberOfTimeSteps:,:,:,:]
x_coriolis = np.copy(myData['tasks']["x_coriolis"])[-numberOfTimeSteps:,:,:,:]
x_inertia = np.copy(myData['tasks']["x_inertia"])[-numberOfTimeSteps:,:,:,:]

# =========================================================================
# Momentum Y equation
# =========================================================================

y_pressure = np.copy(myData['tasks']["y_pressure"])[-numberOfTimeSteps:,:,:,:]
y_diffusion = np.copy(myData['tasks']["y_diffusion"])[-numberOfTimeSteps:,:,:,:]
y_coriolis = np.copy(myData['tasks']["y_coriolis"])[-numberOfTimeSteps:,:,:,:]
y_inertia =	np.copy(myData['tasks']["y_inertia"])[-numberOfTimeSteps:,:,:,:]

# =========================================================================
# Momentum Z equation
# =========================================================================

z_pressure = np.copy(myData['tasks']["z_pressure"])[-numberOfTimeSteps:,:,:,:]
z_diffusion = np.copy(myData['tasks']["z_diffusion"])[-numberOfTimeSteps:,:,:,:]
z_inertia = np.copy(myData['tasks']["z_inertia"])[-numberOfTimeSteps:,:,:,:]
z_buoyancy = np.copy(myData['tasks']["z_bouyancy"])[-numberOfTimeSteps:,:,:,:]

# =========================================================================
# Vorticity X equation
# =========================================================================

vorticity_x_diffusion = np.copy(myData['tasks']["vorticity_x_diffusion"])[-numberOfTimeSteps:,:,:,:]
vorticity_x_coriolis = np.copy(myData['tasks']["vorticity_x_coriolis"])[-numberOfTimeSteps:,:,:,:]
vorticity_x_bouyancy = np.copy(myData['tasks']["vorticity_x_bouyancy"])[-numberOfTimeSteps:,:,:,:]
vorticity_x_inertia = np.copy(myData['tasks']["vorticity_x_inertia"])[-numberOfTimeSteps:,:,:,:]

# =========================================================================
# Vorticity Y equation
# =========================================================================

vorticity_y_diffusion = np.copy(myData['tasks']["vorticity_y_diffusion"])[-numberOfTimeSteps:,:,:,:]
vorticity_y_coriolis = np.copy(myData['tasks']["vorticity_y_coriolis"])[-numberOfTimeSteps:,:,:,:]
vorticity_y_bouyancy = np.copy(myData['tasks']["vorticity_y_bouyancy"])[-numberOfTimeSteps:,:,:,:]
vorticity_y_inertia =	np.copy(myData['tasks']["vorticity_y_inertia"])[-numberOfTimeSteps:,:,:,:]

# =========================================================================
# Vorticity Z equation
# =========================================================================

vorticity_z_diffusion = np.copy(myData['tasks']["vorticity_z_diffusion"])[-numberOfTimeSteps:,:,:,:]
vorticity_z_inertia = np.copy(myData['tasks']["vorticity_z_inertia"])[-numberOfTimeSteps:,:,:,:]
vorticity_z_coriolis = np.copy(myData['tasks']["vorticity_z_coriolis"])[-numberOfTimeSteps:,:,:,:]

myData.close()


hf = h5py.myData('{}{}/last-{}-time-steps.h5'.format(dataDirectory,snapshot_myData_name, str(numberOfTimeSteps)), 'w')

hf.create_dataset('u', data = u)
hf.create_dataset('v', data = v)
hf.create_dataset('w', data = w)
hf.create_dataset('T', data = T)
hf.create_dataset('p', data = p)
hf.create_dataset('uz', data = uz)
hf.create_dataset('vz', data = vz)
hf.create_dataset('wz', data = wz)
hf.create_dataset('Tz', data = Tz)

hf.create_dataset('x_pressure', data = x_pressure)
hf.create_dataset('x_diffusion', data = x_diffusion)
hf.create_dataset('x_coriolis', data = x_coriolis)
hf.create_dataset('x_inertia', data = x_inertia)

hf.create_dataset('y_pressure', data = y_pressure)
hf.create_dataset('y_diffusion', data = y_diffusion)
hf.create_dataset('y_coriolis', data = y_coriolis)
hf.create_dataset('y_inertia', data = y_inertia)

hf.create_dataset('z_pressure', data = z_pressure)
hf.create_dataset('z_diffusion', data = z_diffusion)
hf.create_dataset('z_buoyancy', data = z_buoyancy)
hf.create_dataset('z_inertia', data = z_inertia)

hf.create_dataset('vorticity_x_diffusion', data = vorticity_x_diffusion)
hf.create_dataset('vorticity_x_coriolis', data = vorticity_x_coriolis)
hf.create_dataset('vorticity_x_bouyancy', data = vorticity_x_bouyancy)
hf.create_dataset('vorticity_x_inertia', data = vorticity_x_inertia)

hf.create_dataset('vorticity_y_diffusion', data = vorticity_y_diffusion)
hf.create_dataset('vorticity_y_coriolis', data = vorticity_y_coriolis)
hf.create_dataset('vorticity_y_bouyancy', data = vorticity_y_bouyancy)
hf.create_dataset('vorticity_y_inertia', data = vorticity_y_inertia)

hf.create_dataset('vorticity_z_diffusion', data = vorticity_z_diffusion)
hf.create_dataset('vorticity_z_coriolis', data = vorticity_z_coriolis)
hf.create_dataset('vorticity_z_inertia', data = vorticity_z_inertia)

hf.close()
























