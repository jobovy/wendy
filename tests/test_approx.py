# test_approx.py: some basic tests of the approximate N-body code
import numpy
import pytest
import wendy
numpy.random.seed(2)
def test_energy_conservation():
    # Test that energy is conserved for a simple problem
    x= numpy.array([-1.1,0.1,1.3])
    v= numpy.array([3.,2.,-5.])
    m= numpy.array([1.,1.,1.])
    g= wendy.nbody(x,v,m,0.05,approx=True,nleap=100000)
    E= wendy.energy(x,v,m)
    cnt= 0
    while cnt < 100:
        tx,tv= next(g)
        assert numpy.fabs(wendy.energy(tx,tv,m)-E)/E < 10.**-6., "Energy not conserved during approximate N-body integration"
        cnt+= 1
    return None

def test_energy_conservation_unequalmasses():
    # Test that energy is conserved for a simple problem
    x= numpy.array([-1.1,0.1,1.3])
    v= numpy.array([3.,2.,-5.])
    m= numpy.array([1.,2.,3.])
    g= wendy.nbody(x,v,m,0.05,approx=True,nleap=100000)
    E= wendy.energy(x,v,m)
    cnt= 0
    while cnt < 100:
        tx,tv= next(g)
        assert numpy.fabs(wendy.energy(tx,tv,m)-E)/E < 10.**-6., "Energy not conserved during approximate N-body integration"
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
    g= wendy.nbody(x,v,m,0.05,approx=True,nleap=1000)
    E= wendy.energy(x,v,m)
    cnt= 0
    while cnt < 100:
        tx,tv= next(g)
        assert numpy.fabs(wendy.energy(tx,tv,m)-E)/E < 10.**-6., "Energy not conserved during approximate N-body integration"
        cnt+= 1
    return None

def test_momentum_conservation_unequalmasses():
    # Test that momentum is conserved for a simple problem
    x= numpy.array([-1.1,0.1,1.3])
    v= numpy.array([3.,2.,-5.])
    m= numpy.array([1.,2.,3.])
    g= wendy.nbody(x,v,m,0.05,approx=True,nleap=1000)
    p= wendy.momentum(v,m)
    cnt= 0
    while cnt < 100:
        tx,tv= next(g)
        assert numpy.fabs(wendy.momentum(tv,m)-p) < 10.**-10., "Momentum not conserved during approximate N-body integration"
        cnt+= 1
    return None

def test_notracermasses():
    # approx should work with tracer sheets
    # Test that energy is conserved for a self-gravitating disk
    N= 101
    totmass= 1.
    sigma= 1.
    zh= 2.*sigma**2./totmass
    x= numpy.arctanh(2.*numpy.random.uniform(size=N)-1)*zh
    v= numpy.random.normal(size=N)*sigma
    v-= numpy.mean(v) # stabilize
    m= numpy.ones_like(x)/N*(1.+0.1*(2.*numpy.random.uniform(size=N)-1))
    m[N//2:]= 0.
    m*= 2.
    g= wendy.nbody(x,v,m,0.05,approx=True,nleap=1000)
    E= wendy.energy(x,v,m)
    cnt= 0
    while cnt < 100:
        tx,tv= next(g)
        assert numpy.fabs(wendy.energy(tx,tv,m)-E)/E < 10.**-6., "Energy not conserved during approximate N-body integration with some tracer particles"
        cnt+= 1
    return None

def test_nleap_error():
    # Code should raise ValueError if nleap is not specified for approx. calc.
    x= numpy.array([-1.,1.])
    v= numpy.array([0.,0.])
    m= numpy.array([1.,1.])
    g= wendy.nbody(x,v,m,2,approx=True)
    with pytest.raises(ValueError) as excinfo:
        tx,tv,ncoll, _= next(g)
    assert str(excinfo.value) == 'When approx is True, the number of leapfrog steps nleap= per output time step needs to be set'   
    return None

def test_time():
    # Just run the timer...
    N= 101
    totmass= 1.
    sigma= 1.
    zh= 2.*sigma**2./totmass
    x= numpy.arctanh(2.*numpy.random.uniform(size=N)-1)*zh
    v= numpy.random.normal(size=N)*sigma
    v-= numpy.mean(v) # stabilize
    m= numpy.ones_like(x)/N*(1.+0.1*(2.*numpy.random.uniform(size=N)-1))
    g= wendy.nbody(x,v,m,0.05,approx=True,nleap=1000,full_output=True)
    tx,tv, time_elapsed= next(g)
    assert time_elapsed < 1., 'More than 1 second elapsed for simple problem'
    return None
