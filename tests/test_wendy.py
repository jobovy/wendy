# test_wendy.py: some basic tests of the exact N-body code
import numpy
import pytest
import wendy
numpy.random.seed(2)
def test_energy_conservation():
    # Test that energy is conserved for a simple problem
    x= numpy.array([-1.1,0.1,1.3])
    v= numpy.array([3.,2.,-5.])
    m= numpy.array([1.,1.,1.])
    g= wendy.nbody(x,v,m,0.05)
    E= wendy.energy(x,v,m)
    cnt= 0
    while cnt < 100:
        tx,tv= next(g)
        assert numpy.fabs(wendy.energy(tx,tv,m)-E) < 10.**-10., "Energy not conserved during simple N-body integration"
        cnt+= 1
    return None

def test_energy_conservation_unequalmasses():
    # Test that energy is conserved for a simple problem
    x= numpy.array([-1.1,0.1,1.3])
    v= numpy.array([3.,2.,-5.])
    m= numpy.array([1.,2.,3.])
    g= wendy.nbody(x,v,m,0.05)
    E= wendy.energy(x,v,m)
    cnt= 0
    while cnt < 100:
        tx,tv= next(g)
        assert numpy.fabs(wendy.energy(tx,tv,m)-E) < 10.**-10., "Energy not conserved during simple N-body integration"
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
    g= wendy.nbody(x,v,m,0.05)
    E= wendy.energy(x,v,m)
    cnt= 0
    while cnt < 100:
        tx,tv= next(g)
        assert numpy.fabs(wendy.energy(tx,tv,m)-E) < 10.**-10., "Energy not conserved during simple N-body integration"
        cnt+= 1
    return None

def test_energy_conservation_unequalmasses_python():
    # Test that energy is conserved for a simple problem, using Python method
    x= numpy.array([-1.1,0.1,1.3])
    v= numpy.array([3.,2.,-5.])
    m= numpy.array([1.,2.,3.])
    g= wendy.nbody_python(x,v,m,0.05)
    E= wendy.energy(x,v,m)
    cnt= 0
    while cnt < 100:
        tx,tv= next(g)
        assert numpy.fabs(wendy.energy(tx,tv,m)-E) < 10.**-10., "Energy not conserved during simple N-body integration"
        cnt+= 1
    return None

def test_momentum_conservation_unequalmasses():
    # Test that momentum is conserved for a simple problem
    x= numpy.array([-1.1,0.1,1.3])
    v= numpy.array([3.,2.,-5.])
    m= numpy.array([1.,2.,3.])
    g= wendy.nbody(x,v,m,0.05)
    p= wendy.momentum(v,m)
    cnt= 0
    while cnt < 100:
        tx,tv= next(g)
        assert numpy.fabs(wendy.momentum(tv,m)-p) < 10.**-10., "Momentum not conserved during simple N-body integration"
        cnt+= 1
    return None

def test_energy_individual():
    # Simple test that the individual energies are calculated correctly
    x= numpy.array([-1.1,0.1,1.3])
    v= numpy.array([3.,2.,-5.])
    m= numpy.array([1.,2.,3.])
    E= wendy.energy(x,v,m,individual=True)
    assert numpy.fabs(E[0]-m[0]*v[0]**2./2.-m[0]*(m[1]*numpy.fabs(x[0]-x[1])
                                                  +m[2]*numpy.fabs(x[0]-x[2]))) < 10.**-10
    assert numpy.fabs(E[1]-m[1]*v[1]**2./2.-m[1]*(m[0]*numpy.fabs(x[0]-x[1])
                                                  +m[2]*numpy.fabs(x[2]-x[1]))) < 10.**-10
    assert numpy.fabs(E[2]-m[2]*v[2]**2./2.-m[2]*(m[0]*numpy.fabs(x[0]-x[2])
                                                  +m[1]*numpy.fabs(x[2]-x[1]))) < 10.**-10
    return None

def test_potential():
    # Simple test that the potential is calculated correctly
    x= numpy.array([-1.1,0.1,1.3])
    v= numpy.array([3.,2.,-5.])
    m= numpy.array([1.,2.,3.])
    y= numpy.arange(-2.,2.5,0.5)
    p= wendy.potential(y,x,v,m)
    for ty,tp in zip(y,p):
        assert numpy.fabs(tp-numpy.sum(m*numpy.fabs(x-ty))) < 10.**-10., 'Potential is computed incorrectly'
    return None
                                                                  
def test_count_ncoll():
    # Simple test where we know the number of collisions
    x= numpy.array([-1.,1.])
    v= numpy.array([0.,0.])
    m= numpy.array([1.,1.]) # First collision at t=sqrt(2)
    g= wendy.nbody(x,v,m,1,full_output=True)
    tx,tv,ncoll, _= next(g)
    assert ncoll == 0, 'Number of collisions in simple 2-body problem is wrong'
    tx,tv,ncoll, _= next(g) # collision should have happened now
    assert ncoll == 1, 'Number of collisions in simple 2-body problem is wrong'
    # Next collision is at dt = 2sqrt(2) => dt=2.82 ==> t =~ 4.24
    tx,tv,ncoll, _= next(g)
    assert ncoll == 1, 'Number of collisions in simple 2-body problem is wrong'
    tx,tv,ncoll, _= next(g)
    assert ncoll == 1, 'Number of collisions in simple 2-body problem is wrong'
    tx,tv,ncoll, _= next(g) # collision should have happened now
    assert ncoll == 2, 'Number of collisions in simple 2-body problem is wrong'
    return None

def test_maxncoll_error():
    # Simple test where we know the number of collisions
    x= numpy.array([-1.,1.])
    v= numpy.array([0.,0.])
    m= numpy.array([1.,1.]) # First collision at t=sqrt(2)
    g= wendy.nbody(x,v,m,2,maxcoll=0,full_output=True)
    with pytest.raises(RuntimeError) as excinfo:
        tx,tv,ncoll, _= next(g)
    assert str(excinfo.value) == 'Maximum number of collisions per time step exceeded'   
    return None

def test_maxncoll_warn():
    # Simple test where we know the number of collisions
    x= numpy.array([-1.,1.])
    v= numpy.array([0.,0.])
    m= numpy.array([1.,1.]) # First collision at t=sqrt(2)
    g= wendy.nbody(x,v,m,2,maxcoll=0,full_output=True,warn_maxcoll=True)
    with pytest.warns(RuntimeWarning) as record:
        tx,tv,ncoll, _= next(g)
    # check that only one warning was raised
    assert len(record) == 1
    # check that the message matches
    assert record[0].message.args[0] == "Maximum number of collisions per time step exceeded"
    return None

def test_notracermasses_error():
    # Test that an exception is raised when some of the masses are zero 
    # (or just tiny)
    x= numpy.array([-1.,1.])
    v= numpy.array([0.,0.])
    # Test zero
    m= numpy.array([1.,0.])
    g= wendy.nbody(x,v,m,2)
    with pytest.raises(ValueError) as excinfo:
        tx,tv= next(g)
    assert str(excinfo.value) == 'Tiny masses m much smaller than the median mass are not supported, please remove these from the inputs'   
    # Also test tiny
    m= numpy.array([1.,10.**-10.])
    g= wendy.nbody(x,v,m,2)
    with pytest.raises(ValueError) as excinfo:
        tx,tv= next(g)
    assert str(excinfo.value) == 'Tiny masses m much smaller than the median mass are not supported, please remove these from the inputs'   
    return None

