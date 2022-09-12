"""convert_to_vtr: convert Dedalus fields stored in HDF5 files to vtk format for 3D visualization

Usage:
    convert_to_vtr [--fields=<fields> --nt=<nt>] --dir=<dir> [--output_file=<output_file>]

Options:
    --fields=<fields>             Comma separated list of fields to extract from the hdf5 file [default: None]
    --nt=<nt>                     Time index [default: -1]
    --output_file=<output_file>   Output file name
"""
from dedalus.extras import plot_tools
from pathlib import Path
from docopt import docopt
from pyevtk.hl import gridToVTK
import h5py
import numpy as np
import os

H5_FIELD_PATH = 'tasks/'
H5_SCALE_PATH = 'scales/'
H5_DIM_LABEL = 'DIMENSION_LABELS'
H5_STR_DECODE = 'UTF-8'

if __name__ == "__main__":
    args = docopt(__doc__ )
    nt = int(args['--nt'])
    fields = args['--fields']
    dir = str(args['--dir'])

    # =============================================================================
    # Extract the snapshot filename 
    # =============================================================================

    for idx,lines in enumerate(os.listdir(dir)):
        if lines.find('snapshot') != -1:
            snapshotDirectory = os.listdir(dir)[idx]
    
    infile = dir + snapshotDirectory + '/snapshots.h5'
    outfile = dir + dir.split('/')[1]
    print(outfile)
    if fields is None:
        raise ValueError("Must specify fields to copy.")

    fields = fields.split(',')
    print("fields = {}".format(fields))

    print("outfile = {}".format(outfile))

    datafile = h5py.File(infile,"r")

    field_names = [H5_FIELD_PATH+f for f in fields]
    dim_labels = datafile[field_names[0]].attrs[H5_DIM_LABEL][1:]

    if len(dim_labels) != 3:
        raise NotImplementedError("hdf2vtk only supports 3D data.")

    scale_names = []
    print("#####################################")
    for d in dim_labels:
        
        print(H5_SCALE_PATH)
        scale_names.append(H5_SCALE_PATH+d)
    print("#####################################")
    
    #scale_names = [H5_SCALE_PATH+d.decode(H5_STR_DECODE) for d in dim_labels]
    # just get first scale you find...
    grid_scale = list(datafile[scale_names[0]].keys())[0]
    scale_names = [sn+'/'+grid_scale for sn in scale_names]
    x = plot_tools.get_1d_vertices(datafile[scale_names[0]][:])
    y = plot_tools.get_1d_vertices(datafile[scale_names[1]][:])
    z = plot_tools.get_1d_vertices(datafile[scale_names[2]][:])

    cellData = {}
    for i, f in enumerate(fields):
        #cellData[f] = np.asfortranarray(datafile[field_names[i]][nt])
        cellData[f] = datafile[field_names[i]][nt]

    print(outfile)

    gridToVTK(outfile, x, y, z, cellData = cellData)
