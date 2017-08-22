# test_massless.py: test the treatment of test particles
import numpy
import wendy
numpy.random.seed(2)
def test_basic():
    # Basic test where I've computed the correct evolution for a short time
    x= numpy.array([1.,-1.])
    v= numpy.array([0.,0.])
    m= numpy.array([1.,1.]) # First collision at t=sqrt(2)
    xp= numpy.array([-3.,5.,6.])
    vp= numpy.array([0.,0.,0.])
    g= wendy.nbody(x,v,m,numpy.sqrt(2.)-0.0001,full_output=True,xt=xp,vt=vp)
    # First collision: massive--massive at t=sqrt(2)
    tx,tv,txt,tvt,ncoll, _, _= next(g)
    assert (tx[0] > 0.)*(tx[1] < 0.), "Massive--massive collision shouldn't have happened yet"
    # Now go beyond first massive--massive collision
    g= wendy.nbody(tx,tv,m,0.0002,full_output=True,xt=txt,vt=tvt)
    tx,tv,txt,tvt,ncoll, _, _= next(g)
    assert (tx[0] < 0.)*(tx[1] > 0.), "Massive--massive collision should have happened by now"
    # Now go to first massive--massless collision, at 
    tcoll= numpy.amax(numpy.roots([0.5,3.*numpy.sqrt(2.),-1.]))
    g= wendy.nbody(tx,tv,m,tcoll-0.0002,full_output=True,xt=txt,vt=tvt)
    # Should be just before
    tx,tv,txt,tvt,ncoll, _, _= next(g)
    assert (txt[0] < tx[0]), "Massive--massless collision shouldn't have happened yet"
    g= wendy.nbody(tx,tv,m,0.0002,full_output=True,xt=txt,vt=tvt)
    # Should be just after
    tx,tv,txt,tvt,ncoll, _, _= next(g)
    assert (txt[0] > tx[0]), "Massive--massless collision shouldn't have happened yet"
    # Next collision should be massive 1 with massless 2 around t=0.2685
    g= wendy.nbody(tx,tv,m,0.2685-0.01,full_output=True,xt=txt,vt=tvt)
    # Should be just before
    tx,tv,txt,tvt,ncoll, _, _= next(g)
    assert (txt[0] < tx[1]), "Massive--massless collision shouldn't have happened yet"
    g= wendy.nbody(tx,tv,m,0.02,full_output=True,xt=txt,vt=tvt)
    # Should be just after
    tx,tv,txt,tvt,ncoll, _, _= next(g)
    assert (txt[0] > tx[1]), "Massive--massless collision shouldn't have happened yet" 
    return None

#def test_from_pickle():
#    import pickle
#    with open('test.sav','rb') as savefile:
#        x= pickle.load(savefile)
#        v= pickle.load(savefile)
#        m= pickle.load(savefile)
#        xt= pickle.load(savefile)
#        vt= pickle.load(savefile)
#    g= wendy.nbody(x,v,m,0.05,xt=xt,vt=vt)
#    tx,tv, txt, tvt= next(g)
#    return None

def test_selfgravitating():
    # For an equilibrium configuration, the energy of individual test particles
    # should be approximately conserved (not so for a non-equilibrium config)
    # Test that energy is conserved for a self-gravitating disk, by checking
    # that the energy of test particles is conserved at the same level as that
    # of the massive particles
    N= 1001
    totmass= 1.
    sigma= 1.
    zh= 2.*sigma**2./totmass
    x= numpy.arctanh(2.*numpy.random.uniform(size=N)-1)*zh
    v= numpy.random.normal(size=N)*sigma
    v-= numpy.mean(v) # stabilize
    m= numpy.ones_like(x)/N*(1.+0.1*(2.*numpy.random.uniform(size=N)-1))
    # Also generate massless particles
    M= 201
    xt= numpy.arctanh(2.*numpy.random.uniform(size=M)-1)*zh
    vt= numpy.random.normal(size=M)*sigma
    vt-= numpy.mean(vt) # stabilize   
    g= wendy.nbody(x,v,m,0.05,xt=xt,vt=vt)
    E= wendy.energy(x,v,m,individual=True)
    Et= wendy.energy(x,v,m,individual=True,xt=xt,vt=vt) # dummy m=1
    cnt= 0
    while cnt < 100:
        tx,tv, txt, tvt= next(g)
        assert numpy.amax(numpy.fabs((wendy.energy(tx,tv,m,individual=True,xt=txt,vt=tvt)-Et)/Et)) < 2.*numpy.amax(numpy.fabs((wendy.energy(tx,tv,m,individual=True)-E)/E)), "Energy of massless particles not conserved during simple N-body integration of equilibrium disk"
        cnt+= 1
    return None

def test_massless_stuck_left():
    # Test that the case where a massless particle is stuck between two massive
    # particles that are colliding is handled correctly (such a particle 
    # shouldn't be there, but can end up there because of round-off)
    load_arr= numpy.load('test_data/massless_stuck.npz')
    g= wendy.nbody(load_arr['x'],load_arr['v'],load_arr['m'],
                   0.00000197,
                   xt=load_arr['xt'],vt=load_arr['vt'],
                   maxcoll_tp=10000000000000000,
                   maxcoll=10000000000000000)
    E= wendy.energy(load_arr['x'],load_arr['v'],load_arr['m'])
    tx,tv, txt, tvt= next(g)
    assert numpy.fabs(E-wendy.energy(tx,tv,load_arr['m'])) < 10**-10., 'Energy not conserved when integrating N-body system'
    assert numpy.all(True^numpy.isnan(txt)), 'Massless particles fail when integrating N-body system'
    assert numpy.all(True^numpy.isnan(tvt)), 'Massless particles fail when integrating N-body system'
    return None

def test_massless_stuck_right():
    # Same as above, but mirror, such that the particle is stuck in a different manner
    load_arr= numpy.load('test_data/massless_stuck.npz')
    g= wendy.nbody(-load_arr['x'],-load_arr['v'],load_arr['m'],
                    0.00000197,
                    xt=-load_arr['xt'],vt=-load_arr['vt'],
                    maxcoll_tp=10000000000000000,
                    maxcoll=10000000000000000)
    E= wendy.energy(-load_arr['x'],-load_arr['v'],load_arr['m'])
    tx,tv, txt, tvt= next(g)
    assert numpy.fabs(E-wendy.energy(tx,tv,load_arr['m'])) < 10**-10., 'Energy not conserved when integrating N-body system'
    assert numpy.all(True^numpy.isnan(txt)), 'Massless particles fail when integrating N-body system'
    assert numpy.all(True^numpy.isnan(tvt)), 'Massless particles fail when integrating N-body system'
    return None

def test_gravcollapse_manymanyparticles():
    # Run the gravitational collapse example with 100,000 massless particles
    # until just when it collapses, to resolve multiple collisions at the same 
    # time and generally stress test the integration
    N= 1001
    M= 100002
    dx= numpy.pi
    V0= 0.001
    # Sample massive
    x= dx*(numpy.arange(N)-N//2)/N
    v= -V0*numpy.sin(x)
    m= numpy.ones(N)/float(N)
    # Sample massless
    xtp= dx*(numpy.arange(M)-M//2)/M+numpy.random.uniform(size=M)*dx*10.**-9.
    vtp= -V0*numpy.sin(xtp)
    # Setup generator, go to just before collapse, with pretty big dt
    dt= 0.05
    g= wendy.nbody(x,v,m,dt,xt=xtp,vt=vtp,
                   maxcoll=100000000,maxcoll_tp=10000000000,
                   full_output=True)
    E= wendy.energy(x,v,m)
    cnt= 0
    while cnt < 35:
        tx,tv, txt, tvt,ncoll,_,ncoll_tp= next(g)
        tE= wendy.energy(tx,tv,m)
        assert numpy.fabs(E-tE) < 10.**-10., 'Energy of massive particles not conserved during test'
        assert numpy.fabs((numpy.std(tx)-numpy.std(txt))) < 10.**-4., 'Massless particles do not trace massive particles in x'
        assert numpy.fabs((numpy.std(tv)-numpy.std(tvt))) < 10.**-4., 'Massless particles do not trace massive particles in v'
        cnt+=1
    # Setup generator again, take small steps
    dt= 0.001
    g= wendy.nbody(tx,tv,m,dt,xt=txt,vt=tvt,
                   maxcoll=100000000,maxcoll_tp=10000000000,
                   full_output=True)
    E= wendy.energy(x,v,m)
    cnt= 0
    while cnt < 17:
        tx,tv, txt, tvt,ncoll,_,ncoll_tp= next(g)
        tE= wendy.energy(tx,tv,m)
        assert numpy.fabs(E-tE) < 10.**-10., 'Energy of massive particles not conserved during test'
        assert numpy.fabs((numpy.std(tx)-numpy.std(txt))) < 10.**-4., 'Massless particles do not trace massive particles in x'
        assert numpy.fabs((numpy.std(tv)-numpy.std(tvt))) < 10.**-4., 'Massless particles do not trace massive particles in v'
        cnt+=1
    return None
