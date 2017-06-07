# wendy.py: 1D N-body code
import os
import sys
import sysconfig
import warnings
import ctypes
import ctypes.util
import copy
import numpy
from numpy.ctypeslib import ndpointer
#Find and load the library
_lib= None
outerr= None
PY3= sys.version > '3'
if PY3: # pragma: no cover
    _ext_suffix= sysconfig.get_config_var('EXT_SUFFIX')
else:
    _ext_suffix= '.so'
for path in sys.path:
    try:
        _lib = ctypes.CDLL(os.path.join(path,'wendy_c%s' % _ext_suffix))
    except OSError as e:
        if os.path.exists(os.path.join(path,'wendy_c%s' % _ext_suffix)):# pragma: no cover
            outerr= e
        _lib = None
    else:
        break
if _lib is None: #pragma: no cover
    if not outerr is None:
        warnings.warn("wendy_c extension module not loaded, because of error '%s' " % outerr)
    else:
        warnings.warn("wendy_c extension module not loaded, because wendy_c%s image was not found" % _ext_suffix)
    _ext_loaded= False
else:
    _ext_loaded= True
#Set up the C code
ndarrayFlags= ('C_CONTIGUOUS','WRITEABLE')
_wendy_nbody_onestep_func= _lib._wendy_nbody_onestep
_wendy_nbody_onestep_func.argtypes=\
    [ctypes.c_int,
     ndpointer(dtype=numpy.float64,flags=ndarrayFlags),
     ndpointer(dtype=numpy.float64,flags=ndarrayFlags),
     ndpointer(dtype=numpy.float64,flags=ndarrayFlags),
     ndpointer(dtype=numpy.float64,flags=ndarrayFlags),
     ndpointer(dtype=numpy.int32,flags=ndarrayFlags),
     ctypes.POINTER(ctypes.c_int),
     ctypes.POINTER(ctypes.c_double),
     ndpointer(dtype=numpy.float64,flags=ndarrayFlags),
     ctypes.c_double,
     ctypes.c_int,
     ctypes.POINTER(ctypes.c_int),
     ctypes.POINTER(ctypes.c_int),
     ctypes.POINTER(ctypes.c_double)]
_wendy_nbody_harm_onestep_func= _lib._wendy_nbody_harm_onestep
_wendy_nbody_harm_onestep_func.argtypes=\
    [ctypes.c_int,
     ndpointer(dtype=numpy.float64,flags=ndarrayFlags),
     ndpointer(dtype=numpy.float64,flags=ndarrayFlags),
     ndpointer(dtype=numpy.float64,flags=ndarrayFlags),
     ndpointer(dtype=numpy.float64,flags=ndarrayFlags),
     ndpointer(dtype=numpy.int32,flags=ndarrayFlags),
     ctypes.POINTER(ctypes.c_int),
     ctypes.POINTER(ctypes.c_double),
     ndpointer(dtype=numpy.float64,flags=ndarrayFlags),
     ctypes.c_double,
     ctypes.c_int,
     ctypes.POINTER(ctypes.c_int),
     ctypes.POINTER(ctypes.c_int),
     ctypes.POINTER(ctypes.c_double),
     ctypes.c_double]
class array_w_index(ctypes.Structure):
    _fields_= [("idx", ctypes.c_int),
               ("val", ctypes.c_double)]
_wendy_nbody_approx_onestep_func= _lib._wendy_nbody_approx_onestep
_wendy_nbody_approx_onestep_func.argtypes=\
    [ctypes.c_int,
     ctypes.POINTER(array_w_index),
     ndpointer(dtype=numpy.float64,flags=ndarrayFlags),
     ndpointer(dtype=numpy.float64,flags=ndarrayFlags),
     ndpointer(dtype=numpy.float64,flags=ndarrayFlags),
     ndpointer(dtype=numpy.float64,flags=ndarrayFlags),
     ctypes.c_double,
     ctypes.c_int,
     ctypes.c_double,
     ctypes.POINTER(ctypes.c_int),
     ctypes.POINTER(ctypes.c_double),
     ndpointer(dtype=numpy.float64,flags=ndarrayFlags),
     ndpointer(dtype=numpy.float64,flags=ndarrayFlags)]
_wendy_solve_coll_quad_func= _lib._solve_coll_quad
_wendy_solve_coll_quad_func.argtypes=\
    [ctypes.c_double,ctypes.c_double,ctypes.c_double]
_wendy_solve_coll_quad_func.restype= ctypes.c_double
_wendy_solve_coll_harm_func= _lib._solve_coll_harm
_wendy_solve_coll_harm_func.argtypes=\
    [ctypes.c_double,ctypes.c_double,ctypes.c_double,ctypes.c_double]
_wendy_solve_coll_harm_func.restype= ctypes.c_double

class MyQuadPoly:
    """Simple quadratic polynomial class"""
    def __init__(self,coeff):
        self.coeff= coeff
    def solve(self):
        coeff= self.coeff
        mba= -coeff[1]/coeff[2]
        return 0.5*(mba+numpy.sqrt(mba**2.-4.*coeff[0]/coeff[2]))

def nbody(x,v,m,dt,twopiG=1.,omega=None,approx=False,nleap=None,
          maxcoll=100000,warn_maxcoll=False,
          full_output=False):
    """
    NAME:
       nbody
    PURPOSE:
       run an N-body simulation in 1D
    INPUT:
       x - positions [N]
       v - velocities [N]
       m - masses [N]
       dt - time step
       twopiG= (1.) value of 2 \pi G
       omega= (None) if set, frequency of external harmonic oscillator
       approx= (False) if True, solve the dynamics approximately using leapfrog with exact force evaluations
       nleap= (None) when approx == True, number of leapfrog steps for each dt
       maxcoll= (100000) maximum number of collisions to allow in one time step
       warn_maxcoll= (False) if True, do not raise an error when the maximum number of collisions is exceeded, but instead raise a warning and continue after re-arranging the particles
       full_output= (False) if True, also yield diagnostic information: (a) total number of collisions processed up to this iteration (cumulative; only for exact algorithm), (b) time elapsed resolving collisions if approx is False and for integrating the system if approx is True in just this iteration  (*not* cumulative)
    OUTPUT:
       Generator: each iteration returns (x,v) at equally-spaced time intervals
       + diagnostic info if full_output
    HISTORY:
       2017-04-24 - Written - Bovy (UofT/CCA)
       2017-05-23 - Added omega - Bovy (UofT/CCA)
    """
    if approx: # return approximate solver
        if nleap is None:
            raise ValueError('When approx is True, the number of leapfrog steps nleap= per output time step needs to be set')
        for item in _nbody_approx(x,v,m,dt,nleap,omega=omega,
                                  twopiG=twopiG,full_output=full_output):
            yield item
    # Check omega inpput
    if not omega is None and omega*dt > numpy.pi/2.:
        raise ValueError('When omega is set, omega*dt needs to be less than pi/2; please adjust dt')
    # Check that no masses are tiny, because these cause issues
    if numpy.any(m < numpy.median(m)*10.**-8.):
        raise ValueError('Tiny masses m much smaller than the median mass are not supported, please remove these from the inputs')
    # Sort the data in x
    x= copy.copy(x)
    v= copy.copy(v)
    m= twopiG*copy.copy(m)
    if not omega is None:
        m/= omega**2.
    x,v,m,a,sindx,cindx,next_tcoll,tcoll,err= _setup_arrays(x,v,m,omega=omega)
    ncoll_c= ctypes.c_int(0)
    ncoll= 0
    time_elapsed= ctypes.c_double(0)
    # Simulate the dynamics
    while True:
        if omega is None:
            _wendy_nbody_onestep_func(len(x),
                                      x,v,a,m,sindx,ctypes.byref(cindx),
                                      ctypes.byref(next_tcoll),
                                      tcoll,dt,maxcoll,
                                      ctypes.byref(err),ctypes.byref(ncoll_c),
                                      ctypes.byref(time_elapsed))
        else:
            _wendy_nbody_harm_onestep_func(len(x),
                                           x,v,a,m,sindx,ctypes.byref(cindx),
                                           ctypes.byref(next_tcoll),
                                           tcoll,dt,maxcoll,
                                           ctypes.byref(err),
                                           ctypes.byref(ncoll_c),
                                           ctypes.byref(time_elapsed),
                                           omega)
        ncoll+= ncoll_c.value
        if err.value == -2:
            if warn_maxcoll:
                warnings.warn("Maximum number of collisions per time step exceeded",RuntimeWarning)
                # Re-compute the accelerations
                x,v,m,a,sindx,cindx,next_tcoll,tcoll,err=\
                    _setup_arrays(x,v,m,omega=omega)
            else:
                raise RuntimeError("Maximum number of collisions per time step exceeded")
        if full_output:
            yield(x,v,ncoll,time_elapsed.value)
        else:
            yield(x,v)

def _setup_arrays(x,v,m,omega=None):
    sindx= numpy.argsort(x)
    # Keep track of amount of mass above and below and compute acceleration
    mass_below= numpy.cumsum(m[sindx])
    mass_below[-1]= 0.
    mass_below= numpy.roll(mass_below,1)
    mass_above= numpy.cumsum(m[sindx][::-1])[::-1]
    mass_above[0]= 0.
    mass_above= numpy.roll(mass_above,-1)
    a= (mass_above-mass_below)[numpy.argsort(sindx)]
    # Solve for all collisions, using C code for consistency
    tcoll= []
    for xi,vi,ai,xii,vii,aii in zip(x[sindx][:-1],v[sindx][:-1],a[sindx][:-1],
                                    x[sindx][1:],v[sindx][1:],a[sindx][1:]):
        if omega is None:
            tcoll.append(_wendy_solve_coll_quad_func(xi-xii,vi-vii,(ai-aii)/2.))
        else:
            tcoll.append(_wendy_solve_coll_harm_func(xi-xii,vi-vii,ai-aii,omega))
    tcoll= numpy.array(tcoll)
    cindx= ctypes.c_int(numpy.argmin(tcoll))
    next_tcoll= ctypes.c_double(tcoll[cindx])
    # Prepare for C
    err= ctypes.c_int(0)
    #Array requirements
    x= numpy.require(x,dtype=numpy.float64,requirements=['C','W'])
    v= numpy.require(v,dtype=numpy.float64,requirements=['C','W'])
    a= numpy.require(a,dtype=numpy.float64,requirements=['C','W'])
    m= numpy.require(m,dtype=numpy.float64,requirements=['C','W'])
    sindx= numpy.require(sindx,dtype=numpy.int32,requirements=['C','W'])
    tcoll= numpy.require(tcoll,dtype=numpy.float64,requirements=['C','W'])
    return (x,v,m,a,sindx,cindx,next_tcoll,tcoll,err)

def nbody_python(x,v,m,dt,twopiG=1.):
    """
    NAME:
       nbody_python
    PURPOSE:
       run an N-body simulation in 1D, pure Python version
    INPUT:
       x - positions [N]
       v - velocities [N]
       m - masses [N]
       dt - time step
       twopiG= (1.) value of 2 \pi G
    OUTPUT:
       Generator: each iteration returns (x,v) at equally-spaced time intervals
    HISTORY:
       2017-04-24 - Written - Bovy (UofT/CCA)
    """
    # Sort the data in x
    sindx= numpy.argsort(x)
    i= numpy.arange(len(x)) # to un-sort
    x= x[sindx]
    v= v[sindx]
    m= m[sindx]
    i= i[sindx]
    # Keep track of amount of mass above and below and compute acceleration
    mass_below= numpy.cumsum(m)
    mass_below[-1]= 0.
    mass_below= numpy.roll(mass_below,1)
    mass_above= numpy.cumsum(m[::-1])[::-1]
    mass_above[0]= 0.
    mass_above= numpy.roll(mass_above,-1)
    a= twopiG*(mass_above-mass_below)
    # Setup quadratic polynomials to solve for collisions, 
    # only need to save N-1 differences
    poly= []
    for xi,vi,ai,xii,vii,aii in zip(x[:-1],v[:-1],a[:-1],x[1:],v[1:],a[1:]):
        poly.append(MyQuadPoly([xi-xii,vi-vii,(ai-aii)/2.]))
    # Solve for all collisions
    tcoll= numpy.array([p.solve() for p in poly])
    tcoll[(tcoll < 10.**-10.*dt)]= numpy.inf
    cindx= numpy.argmin(tcoll)
    next_tcoll= tcoll[cindx]
    # Simulate the dynamics
    while True:
        t= numpy.zeros_like(x) # rewind to zero
        while next_tcoll < dt:
            # collide, update collided particles
            x[cindx:cindx+2]+=\
                a[cindx:cindx+2]*(next_tcoll-t[cindx:cindx+2])**2./2.\
                +v[cindx:cindx+2]*(next_tcoll-t[cindx:cindx+2])
            v[cindx:cindx+2]+= \
                a[cindx:cindx+2]*(next_tcoll-t[cindx:cindx+2])
            t[cindx:cindx+2]= next_tcoll # update time
            # swap
            i[cindx], i[cindx+1]= i[cindx+1], i[cindx]
            t[cindx], t[cindx+1]= t[cindx+1], t[cindx]
            x[cindx], x[cindx+1]= x[cindx+1], x[cindx]
            v[cindx], v[cindx+1]= v[cindx+1], v[cindx]
            m[cindx], m[cindx+1]= m[cindx+1], m[cindx]
            mass_below[cindx+1]-= m[cindx+1]-m[cindx]
            mass_above[cindx]+= m[cindx+1]-m[cindx]
            # Update accelerations
            a[cindx:cindx+2]=\
                twopiG*(mass_above[cindx:cindx+2]-mass_below[cindx:cindx+2])
            # Update collision times
            uindx= []
            if cindx > 0:
                uindx.append(cindx-1)
                tdt= t[cindx-1]-next_tcoll
                poly[cindx-1]=\
                    MyQuadPoly([x[cindx-1]+a[cindx-1]*tdt**2./2.-v[cindx-1]*tdt-x[cindx],
                                v[cindx-1]-a[cindx-1]*tdt-v[cindx],
                                a[cindx-1]/2.-a[cindx]/2.])
            uindx.append(cindx)
            poly[cindx]=\
                MyQuadPoly([x[cindx]-x[cindx+1],v[cindx]-v[cindx+1],
                            a[cindx]/2.-a[cindx+1]/2.])
            if cindx < len(x)-2:
                uindx.append(cindx+1)
                tdt= t[cindx+2]-next_tcoll
                poly[cindx+1]=\
                    MyQuadPoly([x[cindx+1]-x[cindx+2]-a[cindx+2]*tdt**2./2.+v[cindx+2]*tdt,
                                v[cindx+1]-v[cindx+2]+a[cindx+2]*tdt,
                                a[cindx+1]/2.-a[cindx+2]/2.])
            t_tcoll= numpy.array([poly[ii].solve()+next_tcoll for ii in uindx])
            tcoll[uindx]= t_tcoll
            cindx= numpy.argmin(tcoll)
            next_tcoll= tcoll[cindx]
        # Update all to next snapshot
        x+=a*(dt-t)**2./2.+v*(dt-t)
        v+= a*(dt-t)
        tcoll-= dt
        next_tcoll-= dt
        yindx= numpy.argsort(i)
        yield (x[yindx],v[yindx])

def _nbody_approx(x,v,m,dt,nleap,omega=None,twopiG=1.,full_output=False):
    """
    NAME:
       _nbody_approx
    PURPOSE:
       run an N-body simulation in 1D, using approximate integration w/ exact forces
    INPUT:
       x - positions [N]
       v - velocities [N]
       m - masses [N]
       dt - output time step
       nleap - number of leapfrog steps / output time step
       omega= (None) if set, frequency of external harmonic oscillator
       twopiG= (1.) value of 2 \pi G
       full_output= (False) if True, also yield diagnostic information: (a) time elapsed for integrating this timestep (*not* cumulative)
    OUTPUT:
       Generator: each iteration returns (x,v) at equally-spaced time intervals
    HISTORY:
       2017-06-03 - Written - Bovy (UofT/CCA)
    """
    if omega is None:
        omega2= -1.
    else:
        omega2= omega**2.
    # Prepare for C
    x= copy.copy(x)
    v= copy.copy(v)
    m= twopiG*copy.copy(m)
    cumulmass= numpy.zeros(len(x))
    revcumulmass= numpy.zeros(len(x))
    err= ctypes.c_int(0)
    #Array requirements
    x= numpy.require(x,dtype=numpy.float64,requirements=['C','W'])
    v= numpy.require(v,dtype=numpy.float64,requirements=['C','W'])
    a= numpy.require(numpy.zeros(len(x)),
                     dtype=numpy.float64,requirements=['C','W'])
    m= numpy.require(m,dtype=numpy.float64,requirements=['C','W'])
    cumulmass= numpy.require(cumulmass,
                             dtype=numpy.float64,requirements=['C','W'])
    revcumulmass= numpy.require(revcumulmass,
                                dtype=numpy.float64,requirements=['C','W'])
    xi= (array_w_index * len(x))()
    for ii in range(len(x)):
        xi[ii].idx= ctypes.c_int(ii)
        xi[ii].val= ctypes.c_double(x[ii])
    # Leapfrog integration
    dt_leap= dt/nleap
    time_elapsed= ctypes.c_double(0)
    while True:
        _wendy_nbody_approx_onestep_func(len(x),ctypes.pointer(xi[0]),
                                         x,v,m,a,dt_leap,nleap,omega2,
                                         ctypes.byref(err),
                                         ctypes.byref(time_elapsed),
                                         cumulmass,revcumulmass)
        if full_output:
            yield(x,v,time_elapsed.value)
        else:
            yield (x,v)

def energy(x,v,m,twopiG=1.,individual=False,omega=None):
    """
    NAME:
       energy
    PURPOSE:
       compute the energy of the system or of each particle
    INPUT:
       x - positions [N]
       v - velocities [N]
       m - masses [N]
       twopiG= (1.) value of 2 \pi G
       individual= (False) if True, return each particle's individual energy (note: individual energies don't add up to the system's energy)
       omega= (None) if set, frequency of external harmonic oscillator
    OUTPUT:
       Energy
    HISTORY:
       2017-04-24 - Written - Bovy (UofT/CCA)
       2017-05-10 - Added individual energies - Bovy (UofT/CCA)
    """
    if not omega is None:
        out= m*omega**2.*x**2./2.
    else:
        out= 0.
    if individual:
        return out\
            +twopiG*m\
            *numpy.sum(m*numpy.fabs(x-numpy.atleast_2d(x).T),axis=1)\
            +m*v**2./2.
    else:
        sindx= numpy.argsort(x)
        mass_below= numpy.roll(numpy.cumsum(m[sindx]),1)
        mass_below[0]= 0.
        xmass_below= numpy.roll(numpy.cumsum((m*x)[sindx]),1)
        xmass_below[0]= 0.
        return numpy.sum(out)\
            +twopiG*numpy.sum(m[sindx]*(mass_below*x[sindx]-xmass_below))\
            +numpy.sum(m*v**2./2.)

def momentum(v,m):
    """
    NAME:
       momentum
    PURPOSE:
       compute the momentum
    INPUT:
       v - velocities [N]
       m - masses [N]
    OUTPUT:
       momentum
    HISTORY:
       2017-04-24 - Written - Bovy (UofT/CCA)
    """
    return numpy.sum(m*v)

def potential(y,x,v,m,twopiG=1.,omega=None):
    """
    NAME:
       potential
    PURPOSE:
       compute the gravitational potential at a set of points
    INPUT:
       y - positions at which to compute the potential
       x - positions of N-body particles [N]
       v - velocities of N-body particles [N]
       m - masses of N-body particles [N]
       twopiG= (1.) value of 2 \pi G
       omega= (None) if set, frequency of external harmonic oscillator
    OUTPUT:
       potential(y)
    HISTORY:
       2017-05-12 - Written - Bovy (UofT/CCA)
    """
    if not omega is None:
        out= omega**2.*y**2./2.
    else:
        out= 0.
    return out\
        +twopiG\
        *numpy.sum(m*numpy.fabs(x-numpy.atleast_2d(y).T),axis=1)
