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
        tx,tv= g.next()
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
        tx,tv= g.next()
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
        tx,tv= g.next()
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
        tx,tv= g.next()
        assert numpy.fabs(wendy.momentum(tv,m)-p) < 10.**-10., "Momentum not conserved during simple N-body integration"
        cnt+= 1
    return None
