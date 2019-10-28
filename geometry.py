# A module to perform various geometry operations
import numpy as np
from useful_functions import logger_setup

#---------------------------------------------------------------------------------------
# MINIMUM IMAGE CONVENTION TOOLS
#---------------------------------------------------------------------------------------
def frac_part(x):
    return x - np.floor(x+0.5) # Element-wise floor from numpy

def distance_pbc(v1, v2, U):
    """Distance between N-dim vectors v1 and v2 according to Minimum Imgage defined by matrix U"""
    Uinv = np.linalg.inv(U)
    dif = v1-v2 #joining vector
    trasf = np.dot(Uinv,dif)  #go to normalized cell-units
    unit = frac_part(trasf) #-0.5 to 0.5 interval
    dif=np.dot(U,unit) #minimal vector
    return np.linalg.norm(dif)

def in_cell(v, M):
    """Return True is N-dim v is inside the cell defined by M False otherwise."""
    Minv = np.linalg.inv(M) # Np is row-wise, we want the matrix to be column wise.
    vt = np.dot(Minv, v)
    if  any([xi<0 for xi in vt]) or any([xi>1 for xi in vt]): # Is there a better way to do it?
        return False
    return True

def map2uc(v, U):
    """Map a position in N-dim v back to the fractional position insde the unit cell defined by U"""
    from math import floor
    Uinv = np.linalg.inv(U) 
    vt = np.dot(Uinv, v)  # Go to normalized cell-units
    vt_unit = np.array([ xi - floor(xi) for xi in vt])  
    v_unit = np.dot(U, vt_unit)
    return v_unit

#---------------------------------------------------------------------------------------
# MANIPULATIONS
#---------------------------------------------------------------------------------------
def vector_lin_stretch(vec, stretch, s_dir):
    """Stretch a vector by a given factor. Optionally along a direction.

    As input: the vector to strerch, the scaling factor, the direction of scaling.
    If direction is not give, the vector is multiplied by the factor (i.e. np_vect*factor)
    """
    c_log = logger_setup(__name__)
    
    norm = np.linalg.norm(s_dir)
    vec = np.array(vec)
    
    # Stretch of a in s dir: v_{stretched}=v_{orth-s}+a*(v.s/|s|)s/|s|
    if norm:
        c_log.info("Directional case: stretching along given direction ("+"%.5f"*len(s)+")", *s)
        s_dir = np.array(s_dir)/norm   # normalize direction, safer
        return vec + (stretch-1.)*dot(vec, s_dir)*s_dir

def vector_iso_stretch(vec, stretch):
    # Isotropic direction
    c_log.info("Isotropic case: stretch along all directions of factor %.5f", stretch)
    return vec*stretch

def pbc_displ(lattice, v):
    """Traslate a layer of a given displacment vector and map it back to unit cell.
    Lattice must be ASE object"""

    from ase import Atoms
    # Check it makes sense
    if type(lattice) != Atoms:
        raise TypeError
    if type(v) != np.array:
        try:
            v = np.array(v)
        except:
            raise ValueError
    
    d_lat = lattice.copy()
    # Move each lattice point. Map it back to cell if needed.
    d_lat.set_positions([map2uc(atm.position+v, lattice.cell.T)
                         for atm in lattice])
    return d_lat

def zcut_geom(geom, z_cut):
    """Split an ASE Object along a plane perpendicular to the vertical (last) dimension"""
    
    from ase import Atoms
    import numpy as np
    from numpy import r_, c_
    c_log = logger_setup(__name__)
    
    # Split top and bottom
    # It's a new cell, starting at origin, positions must be reset
    bottom = Atoms([atm for atm in geom if atm.position[-1] < z_cut])
    top = Atoms([atm for atm in geom if atm.position[-1] > z_cut])

    # Check that it makes sense: no empty objects
    if len(top) == 0 or len(bottom) == 0:
        c_log.error("Cutting plane defines an empty obj: n_top=%i, n_bottom=%i", len(top), len(bottom))
        raise ValueError
    
    # Set top layer positions and cell
    top.set_positions([x-r_[geom.cell[-1][:-1], z_cut]
                      for x in top.positions])
    top_cell = geom.get_cell()
    top_cell[-1] = geom.cell[-1] - [0, 0, z_cut]
    top.set_cell(top_cell)
#    top.wrap() # Should not be needed
    
    # Set bottom layer cell. Here positions should be fine already.
    bottom_cell = geom.get_cell()
    bottom_cell = np.array([v for v in [*geom.cell[:-1],
                                        r_[geom.cell[-1][:-1], z_cut]]
                           ])
    bottom.set_cell(bottom_cell)
    
    # Return the two objects
    return (top, bottom)

def expand_geom(geom, factor):
    """Expand a ASE Atoms structure isotropically by a given factor. Keeps relative coordinate fixed"""
    geom=geom.copy()
    c_log = logger_setup(__name__)
    if [sum(x) for x in geom.cell]==[0,0,0]:
        c_log.error("Supercell is not defined")
        raise ValueError
    
    tmp = geom.get_scaled_positions()
    geom.cell *= factor
    geom.set_scaled_positions(tmp)

    return geom
### END     

#---------------------------------------------------------------------------------------
# EVALUATE A N-DIM PLANE AT A POINT R
#---------------------------------------------------------------------------------------
def plane_at_r(r, n, p):
    """Return the value of a n-dimensional plane at a given (n-1)-point.

    Takes as input the point of interest r (n-dim array), the normal to the plane (n-dim array) and a point intersect by the plane (n-dim array)
    """
    return -sum([ni*(ri-pi) for ni, ri, pi in zip(n[:-1], r[:-1], p[:-1])])/n[-1]+p[-1]


    

