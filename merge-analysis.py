"""Dedalus simulation of 3d Rayleigh benard rotating convection

Usage:
    merge_analysis.py --dir=<directory>
    merge_analysis.py -h | --help

Options:
    -h --help   Display this help message
    --dir=<directory>  H5 File
"""

import subprocess
from dedalus.tools import post
import pathlib
from docopt import docopt
args=docopt(__doc__)

dir = str(args['--dir'])
post.merge_process_files(dir, cleanup=True)
set_paths = list(pathlib.Path(dir).glob('*.h5'))
post.merge_sets('{}/analysis.h5'.format(dir), set_paths, cleanup=True)
