# test_harm.py: some basic tests of the solution in the presence of a background
import numpy
import wendy
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
                                                                  
