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

myData.close()


hf = h5py.File('{}{}/last-{}-time-steps.h5'.format(dataDirectory,snapshot_file_name, str(numberOfTimeSteps)), 'w')

hf.create_dataset('u', data = u)
hf.create_dataset('v', data = v)
hf.create_dataset('w', data = w)
hf.create_dataset('T', data = T)
hf.create_dataset('p', data = p)
hf.create_dataset('uz', data = uz)
hf.create_dataset('vz', data = vz)
hf.create_dataset('wz', data = wz)
hf.create_dataset('Tz', data = Tz)

hf.close()
























