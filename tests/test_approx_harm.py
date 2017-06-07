# test_harm.py: some basic tests of the solution in the presence of a background for the approximate algorithm
import numpy
import wendy
numpy.random.seed(2)
def test_energy_conservation():
    # Test that energy is conserved for a simple problem
    x= numpy.array([-1.1,0.1,1.3])
    v= numpy.array([3.,2.,-5.])
    m= numpy.array([1.,1.,1.])
    omega= 1.1
    g= wendy.nbody(x,v,m,0.05,omega=omega,approx=True,nleap=100000)
    E= wendy.energy(x,v,m,omega=omega)
    cnt= 0
    while cnt < 100:
        tx,tv= next(g)
        assert numpy.fabs(wendy.energy(tx,tv,m,omega=omega)-E)/E < 10.**-6., "Energy not conserved during approximate N-body integration with external harmonic potential"
        cnt+= 1
    return None

def test_energy_conservation_unequalmasses():
    # Test that energy is conserved for a simple problem
    x= numpy.array([-1.1,0.1,1.3])
    v= numpy.array([3.,2.,-5.])
    m= numpy.array([1.,2.,3.])
    omega= 1.1
    g= wendy.nbody(x,v,m,0.05,omega=omega,approx=True,nleap=100000)
    E= wendy.energy(x,v,m,omega=omega)
    cnt= 0
    while cnt < 100:
        tx,tv= next(g)
        assert numpy.fabs(wendy.energy(tx,tv,m,omega=omega)-E)/E < 10.**-6., "Energy not conserved during approximate N-body integration with external harmonic potential"
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
    g= wendy.nbody(x,v,m,0.05,omega=omega,approx=True,nleap=1000)
    E= wendy.energy(x,v,m,omega=omega)
    cnt= 0
    while cnt < 100:
        tx,tv= next(g)
        assert numpy.fabs(wendy.energy(tx,tv,m,omega=omega)-E)/E < 10.**-6., "Energy not conserved during approximate N-body integration with external harmonic potential"
        cnt+= 1
    return None
