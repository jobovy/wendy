# test_wendy.py: some basic tests
import numpy
import wendy
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
                                                                  

