# test_approx_external.py: some basic tests of the solution in the presence of an external force for the approximate algorithm
import numpy
import wendy
numpy.random.seed(2)
def test_energy_conservation():
    # Test that energy is conserved for a simple problem
    x= numpy.array([-1.1,0.1,1.3])
    v= numpy.array([3.,2.,-5.])
    m= numpy.array([1.,1.,1.])
    # Use harmonic oscillator
    omega= 1.1
    ext_force= lambda x,t: -omega**2.*x
    g= wendy.nbody(x,v,m,0.05,ext_force=ext_force,approx=True,nleap=100000)
    E= wendy.energy(x,v,m,omega=omega)
    cnt= 0
    while cnt < 100:
        tx,tv= next(g)
        assert numpy.fabs(wendy.energy(tx,tv,m,omega=omega)-E)/E < 10.**-6., "Energy not conserved during approximate N-body integration with external harmonic potential"
        cnt+= 1
    return None

def test_energy_conservation_unequalmasses_wforceclass():
    # Test that energy is conserved for a simple problem
    x= numpy.array([-1.1,0.1,1.3])
    v= numpy.array([3.,2.,-5.])
    m= numpy.array([1.,2.,3.])
    omega= 1.1
    # Harmonic oscillator force as a class, which numba can't handle
    class Eforce(object):
        def __init__(self,omega):
            self._omega2= omega**2.
        def __call__(self,x,t):
            return -self._omega2*x
    ext_force= Eforce(omega)
    g= wendy.nbody(x,v,m,0.05,ext_force=ext_force,approx=True,nleap=10000)
    E= wendy.energy(x,v,m,omega=omega)
    cnt= 0
    while cnt < 100:
        tx,tv= next(g)
        assert numpy.fabs(wendy.energy(tx,tv,m,omega=omega)-E)/E < 10.**-5., "Energy not conserved during approximate N-body integration with external harmonic potential"
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
    ext_force= lambda x,t: -omega**2.*x
    g= wendy.nbody(x,v,m,0.05,ext_force=ext_force,approx=True,nleap=1000)
    E= wendy.energy(x,v,m,omega=omega)
    cnt= 0
    while cnt < 100:
        tx,tv= next(g)
        assert numpy.fabs(wendy.energy(tx,tv,m,omega=omega)-E)/E < 10.**-6., "Energy not conserved during approximate N-body integration with external harmonic potential"
        cnt+= 1
    return None

def test_energy_conservation_sech2disk_manyparticles_wforceclass():
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
    # Harmonic oscillator force as a class, which numba can't handle
    class Eforce(object):
        def __init__(self,omega):
            self._omega2= omega**2.
        def __call__(self,x,t):
            return -self._omega2*x
    ext_force= Eforce(omega)
    g= wendy.nbody(x,v,m,0.05,ext_force=ext_force,approx=True,nleap=1000)
    E= wendy.energy(x,v,m,omega=omega)
    cnt= 0
    while cnt < 100:
        tx,tv= next(g)
        assert numpy.fabs(wendy.energy(tx,tv,m,omega=omega)-E)/E < 10.**-6., "Energy not conserved during approximate N-body integration with external harmonic potential"
        cnt+= 1
    return None

def test_againstexact_sech2disk_manyparticles():
    # Test that the exact N-body and the approximate N-body agree
    N= 101
    totmass= 1.
    sigma= 1.
    zh= 2.*sigma**2./totmass
    x= numpy.arctanh(2.*numpy.random.uniform(size=N)-1)*zh
    v= numpy.random.normal(size=N)*sigma
    v-= numpy.mean(v) # stabilize
    m= numpy.ones_like(x)/N*(1.+0.1*(2.*numpy.random.uniform(size=N)-1))
    omega= 1.1
    ext_force= lambda x,t: -omega**2.*x
    g= wendy.nbody(x,v,m,0.05,approx=True,nleap=2000,ext_force=ext_force)
    ge= wendy.nbody(x,v,m,0.05,omega=omega)
    cnt= 0
    while cnt < 100:
        tx,tv= next(g)
        txe,tve= next(ge)
        assert numpy.all(numpy.fabs(tx-txe) < 10.**-5.), "Exact and approximate N-body give different positions"
        assert numpy.all(numpy.fabs(tv-tve) < 10.**-5.), "Exact and approximate N-body give different positions"
        cnt+= 1
    return None

