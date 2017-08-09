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
    tx,tv,txt,tvt,ncoll, _= next(g)
    assert (tx[0] > 0.)*(tx[1] < 0.), "Massive--massive collision shouldn't have happened yet"
    # Now go beyond first massive--massive collision
    g= wendy.nbody(tx,tv,m,0.0002,full_output=True,xt=txt,vt=tvt)
    tx,tv,txt,tvt,ncoll, _= next(g)
    assert (tx[0] < 0.)*(tx[1] > 0.), "Massive--massive collision should have happened by now"
    # Now go to first massive--massless collision, at 
    tcoll= numpy.amax(numpy.roots([0.5,3.*numpy.sqrt(2.),-1.]))
    g= wendy.nbody(tx,tv,m,tcoll-0.0002,full_output=True,xt=txt,vt=tvt)
    # Should be just before
    tx,tv,txt,tvt,ncoll, _= next(g)
    assert (txt[0] < tx[0]), "Massive--massless collision shouldn't have happened yet"
    g= wendy.nbody(tx,tv,m,0.0002,full_output=True,xt=txt,vt=tvt)
    # Should be just after
    tx,tv,txt,tvt,ncoll, _= next(g)
    assert (txt[0] > tx[0]), "Massive--massless collision shouldn't have happened yet"
    # Next collision should be massive 1 with massless 2 around t=0.2685
    g= wendy.nbody(tx,tv,m,0.2685-0.01,full_output=True,xt=txt,vt=tvt)
    # Should be just before
    tx,tv,txt,tvt,ncoll, _= next(g)
    assert (txt[0] < tx[1]), "Massive--massless collision shouldn't have happened yet"
    g= wendy.nbody(tx,tv,m,0.02,full_output=True,xt=txt,vt=tvt)
    # Should be just after
    tx,tv,txt,tvt,ncoll, _= next(g)
    assert (txt[0] > tx[1]), "Massive--massless collision shouldn't have happened yet" 
    return None
