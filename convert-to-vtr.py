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

# =============================================================================
# Set some H5 relevant parameters 
# =============================================================================

H5_FIELD_PATH = 'tasks/'
H5_SCALE_PATH = 'scales/'
H5_DIM_LABEL = 'DIMENSION_LABELS'
H5_STR_DECODE = 'UTF-8'

def main():
    """ 
    

    Returns
    -------
    None.
    
    Notes
    -----
    Function which generates a vtr file from h5 data.

    """

    # =============================================================================
    # Extract the docopt arguments 
    # =============================================================================

    args = docopt(__doc__ )
    nt = int(args['--nt'])
    fields = args['--fields']
    directory = str(args['--dir'])    

    # =========================================================================
    # Extract the snapshot filename 
    # =========================================================================

    for idx,lines in enumerate(os.listdir(directory)):
        if lines.find('snapshot') != -1:
            snapshotDirectory = os.listdir(directory)[idx]
    
    # =========================================================================
    # Set the input and output filenames
    # =========================================================================
    
    infile = directory + snapshotDirectory + '/snapshots.h5'
    outfile = directory + directory.split('/')[1]

    if fields is None:
        raise ValueError("Must specify fields to copy.")

    fields = fields.split(',')
    dataFile = h5py.File(infile,"r")

    field_names = [H5_FIELD_PATH+f for f in fields]
    dim_labels = dataFile[field_names[0]].attrs[H5_DIM_LABEL][1:]

    if len(dim_labels) != 3:
        raise NotImplementedError("hdf2vtk only supports 3D data.")

    scale_names = []
    
    print("#"*20)
    for d in dim_labels:
        print(H5_SCALE_PATH)
        scale_names.append(H5_SCALE_PATH+d)
    print("#"*20)

    grid_scale = list(dataFile[scale_names[0]].keys())[0]
    scale_names = [sn+'/'+grid_scale for sn in scale_names]
    x = plot_tools.get_1d_vertices(dataFile[scale_names[0]][:])
    y = plot_tools.get_1d_vertices(dataFile[scale_names[1]][:])
    z = plot_tools.get_1d_vertices(dataFile[scale_names[2]][:])

    cellData = {}
    for i, f in enumerate(fields):
        cellData[f] = dataFile[field_names[i]][nt]


    gridToVTK(outfile, x, y, z, cellData = cellData)


if __name__ == "__main__":
    main()

