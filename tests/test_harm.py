# test_harm.py: some basic tests of the solution in the presence of a background
import numpy
import pytest
import wendy
numpy.random.seed(2)
def test_energy_conservation():
    # Test that energy is conserved for a simple problem
    x= numpy.array([-1.1,0.1,1.3])
    v= numpy.array([3.,2.,-5.])
    m= numpy.array([1.,1.,1.])
    omega= 1.1
    g= wendy.nbody(x,v,m,0.05,omega=omega)
    E= wendy.energy(x,v,m,omega=omega)
    cnt= 0
    while cnt < 100:
        tx,tv= next(g)
        assert numpy.fabs(wendy.energy(tx,tv,m,omega=omega)-E) < 10.**-10., "Energy not conserved during simple N-body integration with external harmonic potential"
        cnt+= 1
    return None

def test_energy_conservation_unequalmasses():
    # Test that energy is conserved for a simple problem
    x= numpy.array([-1.1,0.1,1.3])
    v= numpy.array([3.,2.,-5.])
    m= numpy.array([1.,2.,3.])
    omega= 1.1
    g= wendy.nbody(x,v,m,0.05,omega=omega)
    E= wendy.energy(x,v,m,omega=omega)
    cnt= 0
    while cnt < 100:
        tx,tv= next(g)
        assert numpy.fabs(wendy.energy(tx,tv,m,omega=omega)-E) < 10.**-10., "Energy not conserved during simple N-body integration with external harmonic potential"
        cnt+= 1
    return None

def test_energy_conservation_sech2disk_manyparticles():
    # Test that energy is conserved for a self-gravitating disk
    N= 101
    totmass= 1.
    sigma= 1.
    zh= 2.*sigma**2./totmass
    x= numpy.arctanh(2.*numpy.random.uniform(size=N)-1)*zh
    v= numpy.random.normal(size=N)*sigma
    v-= numpy.mean(v) # stabilize
    m= numpy.ones_like(x)/N*(1.+0.1*(2.*numpy.random.uniform(size=N)-1))
    omega= 1.1
    g= wendy.nbody(x,v,m,0.05,omega=omega)
    E= wendy.energy(x,v,m,omega=omega)
    cnt= 0
    while cnt < 100:
        tx,tv= next(g)
        assert numpy.fabs(wendy.energy(tx,tv,m,omega=omega)-E) < 10.**-10., "Energy not conserved during simple N-body integration with external harmonic potential"
        cnt+= 1
    return None

def test_energy_individual():
    # Simple test that the individual energies are calculated correctly
    x= numpy.array([-1.1,0.1,1.3])
    v= numpy.array([3.,2.,-5.])
    m= numpy.array([1.,2.,3.])
    omega= 1.1
    E= wendy.energy(x,v,m,individual=True,omega=omega)
    assert numpy.fabs(E[0]-m[0]*v[0]**2./2.-m[0]*(m[1]*numpy.fabs(x[0]-x[1])
                                                  +m[2]*numpy.fabs(x[0]-x[2])
                                                  +omega**2.*x[0]**2./2.)) < 10.**-10
    assert numpy.fabs(E[1]-m[1]*v[1]**2./2.-m[1]*(m[0]*numpy.fabs(x[0]-x[1])
                                                  +m[2]*numpy.fabs(x[2]-x[1])
                                                  +omega**2.*x[1]**2./2.)) < 10.**-10
    assert numpy.fabs(E[2]-m[2]*v[2]**2./2.-m[2]*(m[0]*numpy.fabs(x[0]-x[2])
                                                  +m[1]*numpy.fabs(x[2]-x[1])
                                                  +omega**2.*x[2]**2./2.)) < 10.**-10
    return None

def test_potential():
    # Simple test that the potential is calculated correctly
    x= numpy.array([-1.1,0.1,1.3])
    v= numpy.array([3.,2.,-5.])
    m= numpy.array([1.,2.,3.])
    y= numpy.arange(-2.,2.5,0.5)
    omega= 1.1
    p= wendy.potential(y,x,v,m,omega=omega)
    for ty,tp in zip(y,p):
        assert numpy.fabs(tp-numpy.sum(m*numpy.fabs(x-ty))-omega**2.*ty**2./2.) < 10.**-10., 'Potential is computed incorrectly'
    return None
                                                                  
def test_omegadt_lt_pi2():
    # Test that wendy.nbody raises a ValueError when omega x dt > pi/2
    omega= 1.+10.**-8.
    dt= numpy.pi/2.
    with pytest.raises(ValueError) as excinfo:
        g= wendy.nbody(None,None,None,dt,omega=omega)
        next(g)
    assert str(excinfo.value) == 'When omega is set, omega*dt needs to be less than pi/2; please adjust dt'
    return None

def test_maxncoll_error():
    # Simple test where we know the number of collisions
    x= numpy.array([-1.,1.])
    v= numpy.array([0.,0.])
    m= numpy.array([1.,1.]) # First collision at t=~sqrt(2)
    omega= 10.**-2. # to not change the dynamics too much
    g= wendy.nbody(x,v,m,2,maxcoll=0,full_output=True,omega=omega)
    with pytest.raises(RuntimeError) as excinfo:
        tx,tv,ncoll, _= next(g)
    assert str(excinfo.value) == 'Maximum number of collisions per time step exceeded'   
    return None

def test_maxncoll_warn():
    # Simple test where we know the number of collisions
    x= numpy.array([-1.,1.])
    v= numpy.array([0.,0.])
    m= numpy.array([1.,1.]) # First collision at t=~sqrt(2)
    omega= 10.**-2. # to not change the dynamics too much
    g= wendy.nbody(x,v,m,2,maxcoll=0,full_output=True,warn_maxcoll=True,
                   omega=omega)
    with pytest.warns(RuntimeWarning) as record:
        tx,tv,ncoll, _= next(g)
    # check that only one warning was raised
    assert len(record) == 1
    # check that the message matches
    assert record[0].message.args[0] == "Maximum number of collisions per time step exceeded"
    return None

