#!/usr/bin/env python3

"""Cut a geometry along a plane"""

################################################################################
# Preliminaries
################################################################################
import sys, argparse, logging
import numpy as np
from numpy import c_, r_
import ase.io
from ase import Atoms
from useful_functions import logger_setup
from geometry import plane_at_r

def geom_plane_cut(geom, n, p):
    """Return atoms in ASE structure (geom) above and below (tuple of Atoms object) the plane defined by the normal n and intercept p"""

    c_log = logger_setup(__name__)
    
    #-------------------------------------------------------------------------------
    # Divede the points 
    #-------------------------------------------------------------------------------
    p_up, p_down = [], []
    for i, pi in enumerate(geom.positions):
        if pi[-1] > plane_at_r(pi, n, p):
            p_up.append(geom[i])
        else:
            p_down.append(geom[i])
    p_up = Atoms(p_up)
    p_up.set_cell(geom.get_cell())
    p_down = Atoms(p_down)
    p_down.set_cell(geom.get_cell())
    c_log.debug("Up and down")
    c_log.debug(p_up)
    c_log.debug(p_down)

    return p_up, p_down

################################################################################
# Divide atoms with plane 
################################################################################
def plane_cut_wrap(argv):
    """Cut a structure according to given plane

    Valid ASE input geometry from filename or stdin. Plane defined by normal vector n and intercept p.
    Returns a ASE atoms object with the atoms below (or above) the plane.
    If used as script prints an xyz file."""

    #-------------------------------------------------------------------------------
    # Argument parser
    #-------------------------------------------------------------------------------
    parser = argparse.ArgumentParser(description=plane_cut_wrap.__doc__)
    # Positional arguments
    parser.add_argument('filename',
                        default=None,
                        type=str, nargs="?",
                        help='file with initial structure in xzy format. If black use stdin;')
    # Optional args
    parser.add_argument('-n', '--norm',
                        dest='normal',
                        nargs=3, type=float, required=True,
                        help='normal to the plane.')
    parser.add_argument('-p', '--point',
                        dest='point',
                        nargs=3, type=float, required=True,
                        help='intercept of the plane.')
    parser.add_argument('--format',
                        dest='format', default="vasp",
                        help='set ASE-supported format for output.')
    parser.add_argument('-a',
                        action='store_true', dest='get_above',
                        help='get atoms above rather than below the plane.')
    parser.add_argument('--plot',
                        action='store_true', dest='plot_flg',
                        help='plot structure and cutting plane.')
    parser.add_argument('--debug',
                        action='store_true', dest='debug',
                        help='show debug informations.')

    #-------------------------------------------------------------------------------
    # Initialize and check variables
    #-------------------------------------------------------------------------------
    args = parser.parse_args(argv) # Process arguments

    # Set up logger and debug options
    c_log = logger_setup(__name__)
    c_log.setLevel(logging.INFO)
    debug_opt = [] # Pass debug flag down to other functions
    if args.debug:
        c_log.setLevel(logging.DEBUG)
        debug_opt = ['-d']
    c_log.debug(args)

    # Load data from the right source
    # FIXME: there is something broken here. Gets broken pipe.
    if args.filename is None:
        geom = ase.io.read(sys.stdin) #, format="xyz")
    else:
        geom = ase.io.read(args.filename) # , format="xyz")

    # Define the plane
    n = np.array(args.normal)
    n = n/np.linalg.norm(n) # Normalize the normal vector
    p = np.array(args.point)

    #-------------------------------------------------------------------------------
    # Divede the points 
    #-------------------------------------------------------------------------------
    p_up, p_down = geom_plane_cut(geom, n, p)

    xx = np.linspace(min(geom.positions[:,0]), max(geom.positions[:,0]), 5)
    yy = np.linspace(min(geom.positions[:,1]), max(geom.positions[:,1]), 5)
    xm, ym = np.meshgrid(xx, yy)
    zm = -((xm-p[0])*n[0]+(ym-p[1])*n[1])/n[2]+p[2]
    
    #-------------------------------------------------------------------------------
    # Plot the thing
    #-------------------------------------------------------------------------------
    if args.plot_flg:
        import matplotlib.pyplot as plt
        from mpl_toolkits.mplot3d import Axes3D
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        #ax.set_aspect('equal') # Square plot, do not distort angles.
         
        o = np.array([0, 0, 0])
        ax.quiver([0],[0],[0],
                  n[0],n[1],n[2])
        ax.scatter(p[0], p[1], p[2], c="red")


        ax.plot_surface(xm, ym, zm, alpha=0.7)
        
        ax.scatter(p_up.positions[:,0],
                   p_up.positions[:,1],
                   p_up.positions[:,2],
                   alpha=0.8, color="green")
        ax.scatter(p_down.positions[:,0],
                   p_down.positions[:,1],
                   p_down.positions[:,2],
                   alpha=0.8, color="purple")
        
#        ax.set_xlim3d(min(geom.positions[:,0]), max(geom.positions[:,0]))
#        ax.set_xlim3d(min(geom.positions[:,1]), max(geom.positions[:,1]))
#        ax.set_xlim3d(min(geom.positions[:,2]), max(geom.positions[:,2]))
        plt.show()

    # Return points above the plane or below
    if args.get_above:
        res = p_up
    else:
        res = p_down

    if __name__ == "__main__":
        ase.io.write('-', res, format=args.format)

    return res
### End function ---------------------------------------------------------------

################################################################################
# MAIN
################################################################################
if __name__ == "__main__":
    # Bash-script-like functionality: print chemical potentials on stdout and return 0
    plane_cut_wrap(sys.argv[1:])
    exit(0)
